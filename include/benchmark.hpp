// benchmark.hpp - Sistema de benchmarking para Windows/ClangCL
#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <string>
#include <cstdint>
#include <vector>
#include <functional>

namespace bench {

// ============================================================================
// STRUCTS PARA RESULTADOS
// ============================================================================

struct TimingResult {
    double nanoseconds;      // Tiempo en nanosegundos
    uint64_t cpu_cycles;     // Ciclos de CPU (rdtsc)
    double value;           // Valor calculado
    
    // Métodos helper
    double milliseconds() const { return nanoseconds / 1e6; }
    double microseconds() const { return nanoseconds / 1e3; }
    double seconds() const { return nanoseconds / 1e9; }
};

struct BenchmarkStats {
    TimingResult best;       // Mejor ejecución
    TimingResult worst;      // Peor ejecución
    TimingResult average;    // Promedio
    TimingResult median;     // Mediana
    double std_dev;         // Desviación estándar
    int samples;            // Número de muestras
    
    // Métodos para ratios
    double time_ratio(const BenchmarkStats& other) const {
        return average.nanoseconds / other.average.nanoseconds;
    }
    
    double cycles_ratio(const BenchmarkStats& other) const {
        return static_cast<double>(average.cpu_cycles) / other.average.cpu_cycles;
    }
};

// ============================================================================
// CLASE SIMPLIFICADA PARA WINDOWS/CLANGCL
// ============================================================================

class Benchmark {
public:
    // Configuración simple - SIN inicializadores en línea
    struct Config {
        int warmup_iterations;
        int sample_count;
        bool fix_affinity;
        int target_cpu;
        bool verbose;
        
        // Constructor con valores por defecto
        Config() : 
            warmup_iterations(1000),
            sample_count(10),
            fix_affinity(true),
            target_cpu(0),
            verbose(true) {}
    };
    
public:
    // Constructor
    Benchmark(Config config = Config());
    
    // Destructor
    ~Benchmark();
    
    // Método principal
    template<typename Func, typename... Args>
    BenchmarkStats run(const std::string& name, Func&& func, Args&&... args);
    
    // Método de comparación
    template<typename Func1, typename Func2, typename... Args>
    void compare(const std::string& name1, Func1&& func1,
                 const std::string& name2, Func2&& func2,
                 Args&&... args);
    
    // Método estático rápido
    template<typename Func, typename... Args>
    static double quick_time(Func&& func, Args&&... args);
    
    // Utilidades de impresión
    static void print_result(const TimingResult& result, const std::string& name = "");
    static void print_stats(const BenchmarkStats& stats, const std::string& name = "");
    static void print_comparison(const BenchmarkStats& a, const BenchmarkStats& b,
                                const std::string& name_a, const std::string& name_b);
    
private:
    Config config_;
    void* original_affinity_;  // Para restaurar afinidad
    
    // Métodos internos
    void set_cpu_affinity();
    void restore_affinity();
    void warmup();
    uint64_t read_cpu_cycles();
    
    template<typename Func, typename... Args>
    TimingResult measure_single(Func&& func, Args&&... args);
};

} // namespace bench

// Incluir implementación inline
#include "../src/benchmark.inl"  

#endif // BENCHMARK_HPP