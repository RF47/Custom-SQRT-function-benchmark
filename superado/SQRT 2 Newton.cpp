#include <windows.h>
#include <immintrin.h>
#include <limits>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>

using std::cin;
using std::cout;
using namespace std::chrono;


//Versión optimizada con intrinsics para precisión double

// double Sqrt_N(double n) {
//    if (n < 0) {
//        std::cout << "Error: el número no puede ser negativo.\n";
//        return std::numeric_limits<double>::quiet_NaN();  // Mejor manejo de errores
//    }

//    if (n == 0) {
//        return 0;
//    }

//    // Aproximación inicial mejorada
//    double aprox;
//    {
//        uint64_t i;
//        std::memcpy(&i, &n, sizeof(n));
//        i = (i >> 1) + 0x1FF7A30000000000;  // Constante más precisa para double
//        std::memcpy(&aprox, &i, sizeof(i));
//    }

//    // Refinamiento con método de Newton-Raphson optimizado
//    const double error = 1e-15;
//    double prev_aprox;
//    int iterations = 0;
//    const int max_iterations = 6;  // Suficiente para precisión double

//    do {
//        prev_aprox = aprox;
//        aprox = 0.5 * (aprox + n / aprox);
//        iterations++;
//    } while (std::abs(aprox - prev_aprox) > error * aprox && iterations < max_iterations);

//    return aprox;
// }

//static inline double Sqrt_N(double n, int max_iterations = 5) noexcept {
//    if (n <= 0.0) return 0.0;
//    union { double d; uint64_t u; } bits{n};
//    uint64_t exp = (bits.u >> 52) & 0x7FFu;
//    uint64_t half_exp = ((exp - 1023) >> 1) + 1023;
//    bits.u = (bits.u & ~(0x7FFull << 52)) | (half_exp << 52);
//    double ans = bits.d;
//    const double tol = 1e-15;
//    // Cambiamos for(;;) por un for con límite
//    for (int i = 0; i < max_iterations; i++) {
//        double diff = ans * ans - n;
//        if (fabs(diff) <= tol) break;
//        ans = 0.5 * (ans + n / ans);
//    }
//    return ans;
// }
// 

//static inline double Sqrt_N(double n) noexcept {
//    if (n <= 0.0) return 0.0;
//    
//    // Fast path para valores comunes
//    if (n == 1.0) return 1.0;
//    if (n == 4.0) return 2.0;
//    
//    constexpr uint64_t MAGIC = 0x5FE6EB50C7B537A9;
//    uint64_t i;
//    std::memcpy(&i, &n, sizeof(n));
//    i = MAGIC - (i >> 1);
//    
//    double ans;
//    std::memcpy(&ans, &i, sizeof(i));
//    
//    // ans ≈ 1/√n, saltar refinamiento y convertir directamente
//    ans = n * ans;
//    
//    // DOS iteraciones de Newton (suficientes para ~1e-15 error)
//    ans = 0.5 * (ans + n / ans);
//    ans = 0.5 * (ans + n / ans);
//    ans = 0.5 * (ans + n / ans);
//    ans = 0.5 * (ans + n / ans);
//    return ans;
//}

//static inline double Sqrt_N(double n) noexcept {
//    if (n <= 0.0) return 0.0;
//    
//    // Estimación inicial rápida y decente
//    uint64_t i;
//    std::memcpy(&i, &n, sizeof(n));
//    i = (i >> 1) + 0x1FF7A3BEA0000000;
//    
//    double ans;
//    std::memcpy(&ans, &i, sizeof(i));
//    
//    // TRES iteraciones fijas (óptimo balance)
//    ans = 0.5 * (ans + n / ans);  // Iter 1
//    ans = 0.5 * (ans + n / ans);  // Iter 2
//    ans = 0.5 * (ans + n / ans);  // Iter 3
//    
//    return ans;
//}

static inline double Sqrt_N(double n) noexcept {
    if (n <= 0.0) return 0.0;

    // seed para 1/sqrt(n)
    uint64_t i;
    std::memcpy(&i, &n, sizeof(n));
    i = 0x5FE6EB50C7B537A9 - (i >> 1);

    double y;
    std::memcpy(&y, &i, sizeof(i));

    // 2–3 iteraciones suelen alcanzar precisión double
    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);
    y = y * (1.5 - 0.5 * n * y * y);

    return n * y;
}

// // Versión corregida y optimizada con intrinsics
// double Sqrt_N(double n) noexcept {
//    // Manejo de casos especiales
//    if (n < 0.0) {
//        std::cout << "Error: el número no puede ser negativo.\n";
//        return std::numeric_limits<double>::quiet_NaN();
//    }
//    if (n == 0.0) {
//        return 0.0;
//    }
//    // Aproximación inicial usando rsqrt (reciprocal square root) más precisa
//    __m128d x = _mm_set_sd(n);
//    __m128d rsqrt_approx = _mm_set_sd(1.0 / sqrt(n));  // Mejor aproximación inicial
//    // Refinamiento con Newton-Raphson (2 iteraciones suficientes para double)
//    const double half_n = 0.5 * n;
//    double aprox;
//    _mm_store_sd(&aprox, rsqrt_approx);
//    aprox *= n;  // Convertir de 1/sqrt(n) a sqrt(n)
//    // Primera iteración de refinamiento
//    aprox = 0.5 * (aprox + n / aprox);
//    // Segunda iteración (opcional, para máxima precisión)
//    aprox = 0.5 * (aprox + n / aprox);
//    return aprox;
// }

// static double Sqrt_N( const double n) {
//    if (n < 0) {
//        cout << "Error: el número no puede ser negativo.\n";
//        return -1; // Indicar error
//    }
//
//    static double aprox; // Variable para la aproximación
//
//    const double error = 1e-15* n;
//
//    if (n <= 1) {
//        aprox = 1; // Aproximación inicial para números <= 1
//    }
//    else {
//        aprox = 0.5 * n; // Aproximación inicial para números > 1
//    }
//
//     static double n1 = aprox;
//
//    while (abs(n-n1*n1) > error) { // Condición de convergencia
//        n1 = 0.5 * (n1 + n / n1); // Iteración de Newton-Raphson
//
//    }
//    return n1;
//}

//double Sqrt_N(double n) {
//  if (n < 0) {
//    cout << "Error: el número no puede ser negativo.\n";
//    return -1;  // Indicar error
//  }
//
//  if (n == 0) {
//    return 0;
//  }
//
//  // Aproximación inicial utilizando manipulación de bits (segura con memcpy)
//  uint64_t i;
//  double aprox;
//  std::memcpy(&i, &n, sizeof(n));  // Interpretar los bits de `n` como uint64_t
//  i = (i >> 1) + 0x1FF0000000000000;  // Calcular una aproximación inicial
//  std::memcpy(&aprox, &i,
//              sizeof(i));  // Convertir los bits ajustados de nuevo a double
//
//  // Refinamiento con el método de Newton-Raphson
//  const double error = 1e-15;
//  do {
//    aprox = 0.5 * (aprox + n / aprox);
//  } while (abs(aprox * aprox - n) > error * n);
//  return aprox;
//}

int main() {
  SetConsoleOutputCP(CP_UTF8);  // Soporte para UTF-8
  system("color 17");
  cout.precision(15);

  const int r = 100000000;  // Limitar el rango a un valor razonable
  cout << "Test de velocidad hasta " << r << "\n";

  double result = 0.0;

  // Comienza la medición del tiempo para Sqrt_N
  auto start = high_resolution_clock::now();
  for (int i = 1; i <= r; i++) {
    result += Sqrt_N(i);
  }
  auto end = high_resolution_clock::now();
  duration<double, std::micro> elapsed_sqrt_n = end - start;

  cout << "Tiempo de ejecución Sqrt_N: " << elapsed_sqrt_n.count() << " µs\n";
  cout << "Suma Sqrt_N: " << result << "\n";

  // Comienza la medición del tiempo para std::sqrt
  start = high_resolution_clock::now();
  double check = 0.0;
  for (int i = 1; i <= r; i++) {
    check += std::sqrt(i);
  }
  end = high_resolution_clock::now();
  duration<double, std::micro> elapsed_std_sqrt = end - start;
  cout << "Tiempo de ejecución sqrt: " << elapsed_std_sqrt.count() << " µs\n";
  cout << "Suma sqrt: " << check << "\n";

  // Calcular y mostrar el ratio
  double ratio = elapsed_sqrt_n.count() / elapsed_std_sqrt.count();
  cout << "Ratio de tiempos (Sqrt_N / std::sqrt): " << ratio << "\n";

  // Calcular error promedio entre Sqrt_N y std::sqrt
  double error_sum = 0.0;
  for (int i = 1; i <= r; i++) {
    double diff = abs(Sqrt_N(i) - std::sqrt(i));
    error_sum += diff;
  }
  cout << "Error promedio entre Sqrt_N y std::sqrt: " << error_sum / r << "\n";

  // Prueba individual
  double x = 10000000;
  cout << "Raíz cuadrada (Sqrt_N) de " << x << ": " << Sqrt_N(x) << "\n";
  cout << "Raíz cuadrada (sqrt  ) de " << x << ": " << sqrt(x) << "\n";
  cout << "Error absoluto: " << Sqrt_N(x) * Sqrt_N(x) - x << "\n";

  cin.ignore();
  cin.get();
  return 0;
}

