#include "benchmark.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

//  Sqrt_N Custom - Método de Newton-Raphson optimizado
static inline double Sqrt_N(double n) noexcept {
    if (n <= 0.0) return 0.0;

    uint64_t i;
    std::memcpy(&i, &n, sizeof(n));
    i = 0x5FE6EB50C7B537A9 - (i >> 1);

    double y;
    std::memcpy(&y, &i, sizeof(i));

    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);

    return n * y;
}


int main() {

    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    system("chcp 65001 > nul"); 
    // Configurar output
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "========================================\n";
    std::cout << "  BENCHMARK SQRT - Windows/Visual Studio\n";
    std::cout << "========================================\n\n";
    
    // Configuración del benchmark
    bench::Benchmark::Config config;
    config.sample_count = 20;          // Más muestras para mejor estadística
    config.warmup_iterations = 5000;   // Calentamiento adecuado
    config.fix_affinity = true;        // Fijar a un core para consistencia
    config.target_cpu = 0;             // Usar primer core
    config.verbose = true;             // Ver progreso
    
    bench::Benchmark bench(config);
    
    const double test_value = 2.0;
    const double test_values[] = {1.5, 5.0 , 17.0, 101.0, 1000.0, 1000000.0};
    
    // ========================================
    // 1. Benchmark básico
    // ========================================
    std::cout << "\n1. BENCHMARK BÁSICO (valor fijo: " << test_value << ")\n";
    
    auto stats_std = bench.run("std::sqrt", 
        [test_value]() { return std::sqrt(test_value); });
    
    auto stats_custom = bench.run("Sqrt_N (Custom)", 
        [test_value]() { return Sqrt_N(test_value); });
    
    bench::Benchmark::print_comparison(stats_std, stats_custom, 
                                       "std::sqrt", "Sqrt_N");
    
    // ========================================
    // 2. Prueba con diferentes valores (SIMPLE)
    // ========================================
    std::cout << "\n\n2. PRUEBA CON DIFERENTES VALORES\n";
    
    for (double val : test_values) {
        std::cout << "\n--- x = " << std::fixed << std::setprecision(1) << val << " ---\n";
        
        // Calcular y mostrar valores
        double sqrt_std = std::sqrt(val);
        double sqrt_custom = Sqrt_N(val);
        
        std::cout << std::fixed << std::setprecision(15);
        std::cout << "std::sqrt(" << val << ") = " << sqrt_std << "\n";
        std::cout << "Sqrt_N(" << val << ")    = " << sqrt_custom << "\n";
        
        // Mostrar diferencia
        double diff = std::abs(sqrt_custom - sqrt_std);
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Diferencia: " << diff << "\n";
        
        // Medir tiempos
        auto time_std = bench::Benchmark::quick_time([val]() { return std::sqrt(val); });
        auto time_custom = bench::Benchmark::quick_time([val]() { return Sqrt_N(val); });
        
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "Tiempo std::sqrt: " << time_std << " ns\n";
        std::cout << "Tiempo Sqrt_N:    " << time_custom << " ns\n";
        
        if (time_std > 0.1) {
            double ratio = time_custom / time_std;
            std::cout << "Ratio: " << std::setprecision(2) << ratio << "x\n";
        }
    }
    
    // ========================================
    // 3. Benchmark de precisión 
    // ========================================
    std::cout << "\n\n3. PRUEBA DE PRECISIÓN\n";
    
    const int precision_samples = 10000000;
    double max_error = 0.0;
    double avg_error = 0.0;
    double worst_val = 0.0;
    
    for (int i = 1; i <= precision_samples; ++i) {
        double val = 1.0 + (i * 0.01);
        double diff = std::abs(Sqrt_N(val) - std::sqrt(val));
        avg_error += diff;
        
        if (diff > max_error) {
            max_error = diff;
            worst_val = val;
        }
    }
    
    avg_error /= precision_samples;
    
    std::cout << "Muestras: " << precision_samples << "\n";
    
    std::cout << std::scientific << std::setprecision(3); 
    std::cout << "Error máximo:  " << max_error << " (en x=" << worst_val << ")\n";
    std::cout << "Error promedio: " << avg_error << "\n";
    
    // Calcular error relativo también
    double rel_error = max_error / std::sqrt(worst_val);
    std::cout << "Error relativo máximo: " << rel_error << "\n";
    
    // Volver a 2 decimales para el resto
    std::cout << std::fixed << std::setprecision(2);
    
    // ========================================
    // 4. Throughput masivo
    // ========================================
    std::cout << "\n\n4. THROUGHPUT MASIVO\n";
    
    const int massive_iterations = 10000000;
    
    auto start = std::chrono::high_resolution_clock::now();
    volatile double sum_std = 0.0;
    for (int i = 1; i <= massive_iterations; ++i) {
        sum_std += std::sqrt(static_cast<double>(i));
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time_std = std::chrono::duration<double, std::milli>(end - start).count();
    
    start = std::chrono::high_resolution_clock::now();
    volatile double sum_custom = 0.0;
    for (int i = 1; i <= massive_iterations; ++i) {
        sum_custom += Sqrt_N(static_cast<double>(i));
    }
    end = std::chrono::high_resolution_clock::now();
    auto time_custom = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "std::sqrt: " << time_std << " ms, " 
              << (massive_iterations / (time_std / 1000.0) / 1e6) << " M ops/seg\n";
    std::cout << "Sqrt_N:    " << time_custom << " ms, " 
              << (massive_iterations / (time_custom / 1000.0) / 1e6) << " M ops/seg\n";
    std::cout << "Ratio throughput: " << (time_custom / time_std) << "x\n";
    
    // ========================================
    // 5. Información del sistema
    // ========================================
    std::cout << "\n\n5. INFORMACIÓN DEL SISTEMA\n";
    
    SYSTEM_INFO sys_info;
    GetSystemInfo(&sys_info);
    
    std::cout << "Núcleos lógicos: " << sys_info.dwNumberOfProcessors << "\n";
    
    // Obtener frecuencia CPU
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    std::cout << "Frecuencia QPC: " << freq.QuadPart << " Hz\n";
    
    // Medir frecuencia estimada de TSC
    uint64_t tsc_start = __rdtsc();
    Sleep(100);  // 100ms
    uint64_t tsc_end = __rdtsc();
    double tsc_freq = (tsc_end - tsc_start) / 0.1;  // Ciclos por segundo
    std::cout << "Frecuencia TSC estimada: " << (tsc_freq / 1e9) << " GHz\n";
    
    // ========================================
    // 6. Resultado final
    // ========================================
    std::cout << "\n\n========================================\n";
    std::cout << "  RESUMEN FINAL\n";
    std::cout << "========================================\n";
    
    // Calcular throughput de ambos métodos
    double throughput_std = 1e9 / stats_std.average.nanoseconds;  // ops/seg
    double throughput_custom = 1e9 / stats_custom.average.nanoseconds;  // ops/seg
    double throughput_ratio = throughput_custom / throughput_std;
    
    // También usar el throughput masivo como confirmación
    double throughput_massive_std = massive_iterations / (time_std / 1000.0);
    double throughput_massive_custom = massive_iterations / (time_custom / 1000.0);
    double throughput_massive_ratio = throughput_massive_custom / throughput_massive_std;
    
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "THROUGHPUT (benchmark):\n";
    std::cout << "  std::sqrt: " << throughput_std / 1e6 << " M ops/seg\n";
    std::cout << "  Sqrt_N:    " << throughput_custom / 1e6 << " M ops/seg\n";
    
    std::cout << "\nTHROUGHPUT (masivo):\n";
    std::cout << "  std::sqrt: " << throughput_massive_std / 1e6 << " M ops/seg\n";
    std::cout << "  Sqrt_N:    " << throughput_massive_custom / 1e6 << " M ops/seg\n";
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nTIEMPO POR OPERACIÓN:\n";
    std::cout << "  std::sqrt: " << stats_std.average.nanoseconds << " ns\n";
    std::cout << "  Sqrt_N:    " << stats_custom.average.nanoseconds << " ns\n";
    
    // Usar el promedio de ambos ratios para mayor precisión
    double avg_throughput_ratio = (throughput_ratio + throughput_massive_ratio) / 2.0;
    const double SIGNIFICANT_DIFF = 0.05;  // 5% de diferencia significativa
    
    std::cout << "\n📊 CONCLUSIÓN:\n";
    
    if (avg_throughput_ratio > (1.0 + SIGNIFICANT_DIFF)) {
        std::cout << "  ✅ Sqrt_N es " << avg_throughput_ratio 
                  << "x MÁS RÁPIDO que std::sqrt!\n";
    } else if (avg_throughput_ratio < (1.0 - SIGNIFICANT_DIFF)) {
        std::cout << "  ⚠️  Sqrt_N es " << (1.0 / avg_throughput_ratio) 
                  << "x más lento que std::sqrt\n";
    } else {
        std::cout << "  🔄 Rendimiento similar (diferencia < 5%)\n";
        std::cout << "     Ratio: " << avg_throughput_ratio << "x\n";
    }
    
    // Mostrar precisión
    std::cout << "\n🎯 PRECISIÓN:\n";
    std::cout << "  Error promedio: " << std::scientific << std::setprecision(2) 
              << avg_error << "\n";
    std::cout << "  Error máximo:   " << std::scientific << std::setprecision(2) 
              << max_error << "\n";

    std::cout << "\n💡 ANÁLISIS:\n";
    std::cout << "  • Ciclos std::sqrt: " << stats_std.average.cpu_cycles << "\n";
    std::cout << "  • Ciclos Sqrt_N:    " << stats_custom.average.cpu_cycles << "\n";
    int64_t diff = static_cast<int64_t>(stats_custom.average.cpu_cycles) - static_cast<int64_t>(stats_std.average.cpu_cycles);
    std::cout << "  • Diferencia:       " 
              << diff
              << " ciclos\n";
    
    std::cout << "\n========================================\n";
    std::cout << "Presiona Enter para salir...";
    std::cin.ignore();
    std::cin.get();
    
    return 0;
}