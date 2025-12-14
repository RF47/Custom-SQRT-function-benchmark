// benchmark.inl 
#ifndef BENCHMARK_INL
#define BENCHMARK_INL

#include "benchmark.hpp"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <intrin.h>     // Para __rdtsc()
#include <chrono>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <thread>
#include <vector>
#include <iomanip>

// Deshabilitar warnings específicos de ClangCL
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wmissing-braces"

namespace bench {

// ============================================================================
// IMPLEMENTACIÓN Benchmark 
// ============================================================================

inline Benchmark::Benchmark(Config config) 
    : config_(config), original_affinity_(nullptr) {
    
    if (config_.verbose) {
        std::cout << "[Benchmark] Inicializando...\n";
    }
    
    if (config_.fix_affinity) {
        set_cpu_affinity();
    }
    
    warmup();
}

inline Benchmark::~Benchmark() {
    if (config_.fix_affinity) {
        restore_affinity();
    }
}

inline void Benchmark::set_cpu_affinity() {
    HANDLE process = GetCurrentProcess();
    HANDLE thread = GetCurrentThread();
    
    // Guardar afinidad original
    DWORD_PTR process_mask, system_mask;
    GetProcessAffinityMask(process, &process_mask, &system_mask);
    original_affinity_ = reinterpret_cast<void*>(static_cast<uintptr_t>(process_mask));
    
    // Fijar a un solo CPU
    DWORD_PTR new_mask = 1ULL << config_.target_cpu;
    SetThreadAffinityMask(thread, new_mask);
    
    // Alta prioridad
    SetThreadPriority(thread, THREAD_PRIORITY_HIGHEST);
    SetPriorityClass(process, HIGH_PRIORITY_CLASS);
    
    if (config_.verbose) {
        std::cout << "[Benchmark] CPU fijado al core " << config_.target_cpu << "\n";
    }
}

inline void Benchmark::restore_affinity() {
    if (original_affinity_) {
        HANDLE thread = GetCurrentThread();
        DWORD_PTR original = static_cast<DWORD_PTR>(reinterpret_cast<uintptr_t>(original_affinity_));
        SetThreadAffinityMask(thread, original);
        
        // Restaurar prioridad
        SetThreadPriority(thread, THREAD_PRIORITY_NORMAL);
        SetPriorityClass(GetCurrentProcess(), NORMAL_PRIORITY_CLASS);
    }
}

inline void Benchmark::warmup() {
    if (config_.verbose) {
        std::cout << "[Benchmark] Calentando CPU...\n";
    }
    
    volatile double sum = 0.0;
    for (int i = 0; i < config_.warmup_iterations; ++i) {
        sum += sqrt(static_cast<double>(i));
    }
    
    // Pausa más larga para estabilizar
    Sleep(50);
}

inline uint64_t Benchmark::read_cpu_cycles() {
    // Leer Time Stamp Counter con serialización
    _mm_lfence();
    uint64_t cycles = __rdtsc();
    _mm_lfence();
    return cycles;
}

template<typename Func, typename... Args>
TimingResult Benchmark::measure_single(Func&& func, Args&&... args) {
    TimingResult result;
    
    // EJECUTAR MÚLTIPLES VECES para evitar tiempos de 0ns
    const int INNER_ITERATIONS = 100000;  // <-- CLAVE: Ejecutar muchas veces
    
    // Memory barriers
    _ReadWriteBarrier();
    _mm_mfence();
    
    // 1. Medir ciclos
    uint64_t start_cycles = read_cpu_cycles();
    
    // 2. Medir tiempo con QueryPerformanceCounter (más preciso en Windows)
    LARGE_INTEGER start_time, frequency;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start_time);
    
    // 3. Ejecutar función MUCHAS VECES
    volatile double accumulated = 0.0;
    for (int i = 0; i < INNER_ITERATIONS; ++i) {
        accumulated += std::invoke(std::forward<Func>(func), 
                                  std::forward<Args>(args)...);
    }
    
    // 4. Tomar mediciones finales
    LARGE_INTEGER end_time;
    QueryPerformanceCounter(&end_time);
    uint64_t end_cycles = read_cpu_cycles();
    
    // Memory barrier final
    _mm_mfence();
    _ReadWriteBarrier();
    
    // 5. Calcular tiempo promedio por iteración
    double elapsed_seconds = static_cast<double>(end_time.QuadPart - start_time.QuadPart) / 
                            static_cast<double>(frequency.QuadPart);
    
    // Tiempo en NANOSEGUNDOS por iteración
    result.nanoseconds = (elapsed_seconds * 1e9) / INNER_ITERATIONS;
    
    // Ciclos de CPU por iteración
    result.cpu_cycles = (end_cycles - start_cycles) / INNER_ITERATIONS;
    
    // Valor calculado (para evitar optimizaciones)
    result.value = accumulated;
    
    // Filtrar valores imposibles
    if (result.nanoseconds < 0.1) {
        result.nanoseconds = 0.1;  // Mínimo realista
    }
    
    return result;
}

template<typename Func, typename... Args>
BenchmarkStats Benchmark::run(const std::string& name, Func&& func, Args&&... args) {
    if (config_.verbose) {
        std::cout << "[Benchmark] Ejecutando: " << name << "\n";
        std::cout << "[Benchmark] Muestras: " << config_.sample_count << "\n";
    }
    
    std::vector<TimingResult> samples;
    samples.reserve(config_.sample_count);
    
    // Warmup específico (más intensivo)
    for (int i = 0; i < config_.warmup_iterations * 2; ++i) {
        volatile auto dummy = std::invoke(std::forward<Func>(func), 
                                         std::forward<Args>(args)...);
    }
    
    // Tomar muestras
    for (int i = 0; i < config_.sample_count; ++i) {
        if (config_.verbose && (i % 5 == 0 || i == config_.sample_count - 1)) {
            std::cout << "[Benchmark] Muestra " << (i + 1) << "/" 
                      << config_.sample_count << "\n";
        }
        
        samples.push_back(measure_single(std::forward<Func>(func), 
                                        std::forward<Args>(args)...));
        
        // Pausa más larga entre muestras
        if (i < config_.sample_count - 1) {
            Sleep(20);
        }
    }
    
    // ============================================
    // CÁLCULO DE ESTADÍSTICAS MEJORADO
    // ============================================
    
    BenchmarkStats stats;
    stats.samples = config_.sample_count;
    
    if (samples.empty()) {
        return stats;
    }
    
    // 1. Filtrar outliers (valores extremos)
    std::vector<double> valid_times;
    for (const auto& sample : samples) {
        if (sample.nanoseconds > 0.1 && sample.nanoseconds < 1000000.0) {
            valid_times.push_back(sample.nanoseconds);
        }
    }
    
    if (valid_times.empty()) {
        // Si todos son inválidos, usar el primero
        stats.best = samples[0];
        stats.worst = samples[0];
        stats.average = samples[0];
        stats.median = samples[0];
        stats.std_dev = 0.0;
        return stats;
    }
    
    // 2. Ordenar para percentiles
    std::sort(valid_times.begin(), valid_times.end());
    
    // 3. Usar percentil 10 como mejor, percentil 90 como peor (ignora outliers)
    size_t best_idx = valid_times.size() * 10 / 100;
    size_t worst_idx = valid_times.size() * 90 / 100;
    
    if (best_idx >= valid_times.size()) best_idx = 0;
    if (worst_idx >= valid_times.size()) worst_idx = valid_times.size() - 1;
    
    // 4. Encontrar muestras correspondientes a los percentiles
    auto find_closest_sample = [&](double target_ns) -> TimingResult {
        TimingResult closest = samples[0];
        double min_diff = std::abs(samples[0].nanoseconds - target_ns);
        
        for (const auto& sample : samples) {
            double diff = std::abs(sample.nanoseconds - target_ns);
            if (diff < min_diff) {
                min_diff = diff;
                closest = sample;
            }
        }
        return closest;
    };
    
    stats.best = find_closest_sample(valid_times[best_idx]);
    stats.worst = find_closest_sample(valid_times[worst_idx]);
    
    // 5. Calcular promedio usando percentiles 25-75 (ignorar outliers extremos)
    double sum_ns = 0.0;
    uint64_t sum_cycles = 0;
    int count = 0;
    
    size_t q1 = valid_times.size() * 25 / 100;
    size_t q3 = valid_times.size() * 75 / 100;
    
    for (size_t i = q1; i <= q3 && i < valid_times.size(); ++i) {
        for (const auto& sample : samples) {
            if (std::abs(sample.nanoseconds - valid_times[i]) < 0.01) {
                sum_ns += sample.nanoseconds;
                sum_cycles += sample.cpu_cycles;
                count++;
                break;
            }
        }
    }
    
    if (count > 0) {
        stats.average.nanoseconds = sum_ns / count;
        stats.average.cpu_cycles = sum_cycles / count;
    } else {
        // Fallback: promedio simple
        for (const auto& sample : samples) {
            sum_ns += sample.nanoseconds;
            sum_cycles += sample.cpu_cycles;
        }
        stats.average.nanoseconds = sum_ns / samples.size();
        stats.average.cpu_cycles = sum_cycles / samples.size();
    }
    
    // 6. Calcular mediana (percentil 50)
    size_t median_idx = valid_times.size() / 2;
    stats.median = find_closest_sample(valid_times[median_idx]);
    
    // 7. Calcular desviación estándar SOLO de las muestras usadas para el promedio
    double variance = 0.0;
    for (size_t i = q1; i <= q3 && i < valid_times.size(); ++i) {
        double diff = valid_times[i] - stats.average.nanoseconds;
        variance += diff * diff;
    }
    
    if (count > 1) {
        stats.std_dev = std::sqrt(variance / count);
    } else {
        stats.std_dev = 0.0;
    }
    
    return stats;
}

template<typename Func1, typename Func2, typename... Args>
void Benchmark::compare(const std::string& name1, Func1&& func1,
                       const std::string& name2, Func2&& func2,
                       Args&&... args) {
    auto stats1 = run(name1, std::forward<Func1>(func1), std::forward<Args>(args)...);
    auto stats2 = run(name2, std::forward<Func2>(func2), std::forward<Args>(args)...);
    
    print_comparison(stats1, stats2, name1, name2);
}

template<typename Func, typename... Args>
double Benchmark::quick_time(Func&& func, Args&&... args) {
    // Medición rápida sin instanciar Benchmark completo
    // Usar chrono de alta resolución que es más confiable para tiempos cortos
    
    const int ITERATIONS = 100000;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    volatile double sum = 0.0;
    for (int i = 0; i < ITERATIONS; ++i) {
        sum += std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    return static_cast<double>(duration.count()) / ITERATIONS;
}

// ============================================================================
// FUNCIONES DE IMPRESIÓN (MEJORADAS)
// ============================================================================

inline void Benchmark::print_result(const TimingResult& result, const std::string& name) {
    if (!name.empty()) {
        std::cout << "\n" << name << ":\n";
    }
    
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Tiempo:     " << result.nanoseconds << " ns";
    
    if (result.nanoseconds >= 1000.0) {
        std::cout << " (" << result.microseconds() << " µs)";
    }
    std::cout << "\n";
    
    std::cout << "  Ciclos CPU: " << result.cpu_cycles << "\n";
    
    // Calcular frecuencia estimada solo si tiene sentido
    if (result.nanoseconds > 1.0 && result.cpu_cycles > 10) {
        double ghz = (result.cpu_cycles / result.nanoseconds) * 1000.0;
        std::cout << "  CPU estim:  " << std::setprecision(2) << ghz << " GHz\n";
    }
}

inline void Benchmark::print_stats(const BenchmarkStats& stats, const std::string& name) {
    std::cout << "\n════════════════════════════════════════════════════\n";
    if (!name.empty()) {
        std::cout << "  " << name << " (" << stats.samples << " muestras)\n";
    } else {
        std::cout << "  Resultados (" << stats.samples << " muestras)\n";
    }
    std::cout << "════════════════════════════════════════════════════\n";
    
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  Mejor:      " << stats.best.nanoseconds << " ns";
    std::cout << " (" << stats.best.cpu_cycles << " ciclos)\n";
    
    std::cout << "  Peor:       " << stats.worst.nanoseconds << " ns";
    std::cout << " (" << stats.worst.cpu_cycles << " ciclos)\n";
    
    std::cout << "  Promedio:   " << stats.average.nanoseconds << " ns";
    std::cout << " (" << stats.average.cpu_cycles << " ciclos)\n";
    
    std::cout << "  Mediana:    " << stats.median.nanoseconds << " ns";
    std::cout << " (" << stats.median.cpu_cycles << " ciclos)\n";
    
    // Mostrar variación porcentual de forma más clara
    double variation = (stats.std_dev / stats.average.nanoseconds) * 100.0;
    std::cout << "  StdDev:     " << stats.std_dev << " ns";
    std::cout << " (" << std::setprecision(1) << variation << "%)\n";
    
    // Calcular métricas adicionales
    std::cout << "\n  Por iteración:\n";
    std::cout << "    Mejor:    " << stats.best.nanoseconds << " ns/iter\n";
    std::cout << "    Promedio: " << stats.average.nanoseconds << " ns/iter\n";
    std::cout << "    Peor:     " << stats.worst.nanoseconds << " ns/iter\n";
    
    // Throughput estimado (solo si el tiempo es razonable)
    if (stats.average.nanoseconds > 0.1) {
        double throughput = 1e9 / stats.average.nanoseconds;
        std::cout << "\n  Throughput estimado:\n";
        std::cout << "    " << std::setprecision(0) << throughput / 1e6 << " M ops/seg\n";
    }
    
    std::cout << "════════════════════════════════════════════════════\n";
}

inline void Benchmark::print_comparison(const BenchmarkStats& a, const BenchmarkStats& b,
                                       const std::string& name_a, const std::string& name_b) {
    std::cout << "\n════════════════════════════════════════════════════\n";
    std::cout << "  COMPARACIÓN: " << name_a << " vs " << name_b << "\n";
    std::cout << "════════════════════════════════════════════════════\n";
    
    // Calcular throughput de cada función
    double throughput_a = 1e9 / a.average.nanoseconds;  // ops/segundo
    double throughput_b = 1e9 / b.average.nanoseconds;  // ops/segundo
    double throughput_ratio = throughput_b / throughput_a;  // b vs a
    
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  TIEMPO:\n";
    std::cout << "    " << name_a << ": " << a.average.nanoseconds << " ns\n";
    std::cout << "    " << name_b << ": " << b.average.nanoseconds << " ns\n";
    
    // Usar throughput para determinar quién es más rápido
    const double THRESHOLD = 0.02;  // 2% de diferencia
    
    if (throughput_ratio > (1.0 + THRESHOLD)) {
        // b es más rápido que a
        std::cout << "    " << name_b << " es " << std::setprecision(2) 
                  << throughput_ratio << "x más rápido\n";
    } else if (throughput_ratio < (1.0 - THRESHOLD)) {
        // b es más lento que a
        std::cout << "    " << name_b << " es " << std::setprecision(2) 
                  << (1.0 / throughput_ratio) << "x más lento\n";
    } else {
        // Prácticamente iguales
        std::cout << "    " << "Rendimiento similar (diferencia < 2%)\n";
    }
    
    std::cout << "\n  THROUGHPUT:\n";
    std::cout << "    " << name_a << ": " << std::setprecision(0) 
              << throughput_a / 1e6 << " M ops/seg\n";
    std::cout << "    " << name_b << ": " << std::setprecision(0) 
              << throughput_b / 1e6 << " M ops/seg\n";
    
    std::cout << "\n  CICLOS CPU:\n";
    std::cout << "    " << name_a << ": " << a.average.cpu_cycles << " ciclos\n";
    std::cout << "    " << name_b << ": " << b.average.cpu_cycles << " ciclos\n";
    
    double cycles_ratio = static_cast<double>(b.average.cpu_cycles) / a.average.cpu_cycles;
    if (cycles_ratio > (1.0 + THRESHOLD)) {
        std::cout << "    " << name_b << " usa " << std::setprecision(2) 
                  << cycles_ratio << "x más ciclos\n";
    } else if (cycles_ratio < (1.0 - THRESHOLD)) {
        std::cout << "    " << name_b << " usa " << std::setprecision(2) 
                  << (1.0 / cycles_ratio) << "x menos ciclos\n";
    }
    
    std::cout << "\n  ESTABILIDAD:\n";
    double var_a = (a.std_dev / a.average.nanoseconds) * 100.0;
    double var_b = (b.std_dev / b.average.nanoseconds) * 100.0;
    
    std::cout << "    " << name_a << ": " << std::setprecision(1) << var_a << "%\n";
    std::cout << "    " << name_b << ": " << std::setprecision(1) << var_b << "%\n";
    
    std::cout << "════════════════════════════════════════════════════\n";
}

} // namespace bench

#pragma clang diagnostic pop

#endif // BENCHMARK_INL