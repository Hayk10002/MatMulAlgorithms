#include <iostream>
#include <format>
#include <vector>
#include <span>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <functional>
#include <optional>
#include <random>
#include <map>
#include <variant>
#include <format>
#include <ranges>
#include <array>
#include <set>

#include "include/MatrixView.hpp"
#include "include/iterative.hpp"
#include "include/block.hpp"
#include "include/recursive.hpp"
#include "include/MatrixMultiplier.hpp"
#include "include/cmd_args.hpp"

using namespace std::string_view_literals;

int generateRandomNumber(int min = 0, int max = 100) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

template<class F, bool verify>
auto time(F f, int N, int M, int P)
{
    static std::map<std::array<int, 3>, std::array<std::vector<int>, 4>> cache{};
    auto it = cache.find({N, M, P});
    if (it == cache.end())
    {
        std::vector<int> A(N * M);
        std::vector<int> B(M * P);
        std::vector<int> C(N * P);

        std::vector<int> E(N * P);

        std::generate(A.begin(), A.end(), [](){ return generateRandomNumber(); });
        std::generate(B.begin(), B.end(), [](){ return generateRandomNumber(); });
        std::generate(C.begin(), C.end(), [](){ return generateRandomNumber(); });

        MatrixView A_view(A, M);
        MatrixView B_view(B, P);
        MatrixView C_view(C, P);

        MatrixView E_view(E, P);

        if constexpr(verify) 
            naiveCacheFriendlyMatMul(A_view, B_view, E_view, MatMulMode::Add);

        cache[{N, M, P}] = {std::move(A), std::move(B), std::move(C), std::move(E)};
        it = cache.find({N, M, P});
    }

    std::vector<int>& A = it->second[0];
    std::vector<int>& B = it->second[1];
    std::vector<int>& C = it->second[2];
    std::vector<int>& E = it->second[3];

    MatrixView A_view(A, M);
    MatrixView B_view(B, P);
    MatrixView C_view(C, P);

    MatrixView E_view(E, P);

    auto start = std::chrono::high_resolution_clock::now();

    f(A_view, B_view, C_view, MatMulMode::Overwrite);

    auto t_overwrite = std::chrono::high_resolution_clock::now() - start;
    bool overwrite_check = verify ? C == E : true;

    if constexpr(verify) 
        std::ranges::for_each(E, [](int& x){ x *= 2; });

    start = std::chrono::high_resolution_clock::now();

    f(A_view, B_view, C_view, MatMulMode::Add);

    auto t_add = std::chrono::high_resolution_clock::now() - start;
    bool add_check = verify ? C == E : true;

    if constexpr(verify)
        std::ranges::for_each(E, [](int& x){ x /= 2; });

    if constexpr(verify)
        return std::pair{ 
            std::pair{ t_overwrite, overwrite_check }, 
            std::pair{ t_add, add_check} 
        };
    
    else
        return std::pair{ 
            t_overwrite, t_add
        };
}

MatrixMultiplier recursive_until_size(int size)
{
    return MatrixMultiplier::recursive_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
}

MatrixMultiplier strassen_until_size(int size)
{
    return MatrixMultiplier::strassen_then([size](int n, int m, int p){ return n < size || m < size || p < size; }, MatrixMultiplier::naive_cache_friendly_mutliplier);
}

template<class F>
void print_test(std::string_view name, F f, int N, int M, int P, bool verify = true)
{
    static auto to_ms = [](auto t){ return std::chrono::duration_cast<std::chrono::milliseconds>(t); };
    if (verify)
    {
        auto res = time<F, true>(f, N, M, P);
        std::cout << std::format("{}: OWT: {:>7} {:>5}, ADD: {:>7} {:>5}\n", name, to_ms(res.first.first), res.first.second, to_ms(res.second.first), res.second.second);
    }
    else
    {
        auto res = time<F, false>(f, N, M, P);
        std::cout << std::format("{}: OWT: {:>7}, ADD: {:>7}\n", name, to_ms(res.first), to_ms(res.second));
    }
}

struct TestableType
{
    enum class Type
    {
        Naive,
        NaiveMatMul,
        BetterNaive,
        BetterNaiveMatMul,
        Blocked,
        BlockedMatMul,
        Recursive,
        Strassen,
        Hybrid,
        Multithreaded
    } type;
    int val{0};

    using enum Type;

    TestableType(Type type) : type(type) {}
    TestableType(Type type, int val) : type(type), val(val) {} 

    auto operator<=>(const TestableType&) const = default;
    bool operator==(const TestableType&) const = default;
};

struct Testable
{
    using MultiplierType = std::function<void(MatrixView, MatrixView, MatrixView, MatMulMode)>;
    Testable(std::string_view name, auto f) : name(name), f(f) {}
    Testable(TestableType type)
    {
        switch(type.type)
        {
            case TestableType::Naive:             f = naiveMatMul; break;
            case TestableType::NaiveMatMul:       f = naiveMatMul; break;
            case TestableType::BetterNaive:       f = naiveCacheFriendlyMatMul; break;
            case TestableType::BetterNaiveMatMul: f = naiveCacheFriendlyMatMul; break;
            case TestableType::Blocked:           f = cacheFriendlyBlockMatMul; break;
            case TestableType::BlockedMatMul:     f = cacheFriendlyBlockMatMul; break;
            case TestableType::Recursive:         f = recursive_until_size(type.val); break;
            case TestableType::Strassen:          f = strassen_until_size(type.val); break;
            case TestableType::Hybrid:            f = MatrixMultiplier::hybrid_multiplier; break;
            case TestableType::Multithreaded:     f = MatrixMultiplier::multithreaded_hybrid_multiplier; break;
        }

        name = name_from_type(type);
    }

    static std::string short_name_from_type(TestableType type)
    {
        switch(type.type)
        {
            case TestableType::Naive:             return "naive";
            case TestableType::NaiveMatMul:       return "naive_MatMul";
            case TestableType::BetterNaive:       return "better_naive";
            case TestableType::BetterNaiveMatMul: return "better_naive_MatMul";
            case TestableType::Blocked:           return "blocked";
            case TestableType::BlockedMatMul:     return "blocked_MatMul";
            case TestableType::Recursive:         return std::vformat("recursive{}", std::make_format_args(type.val)); break;
            case TestableType::Strassen:          return std::vformat("strassen{}", std::make_format_args(type.val)); break;
            case TestableType::Hybrid:            return "hybrid";
            case TestableType::Multithreaded:     return "multithreaded";
        }
        return "";
    }

    template<class T = int>
    static std::string name_from_type(TestableType type, std::optional<T> placeholder = std::nullopt)
    {
        std::string val = std::format("{:<3}", type.val);
        if (placeholder && type.val == 0) val = std::format("{}", *placeholder);
        switch(type.type)
        {
            case TestableType::Naive:             return             "Naive                                  ";
            case TestableType::NaiveMatMul:       return             "Naive                (MatrixMultiplier)";
            case TestableType::BetterNaive:       return             "Cache-friendly naive                   ";
            case TestableType::BetterNaiveMatMul: return             "Cache-friendly naive (MatrixMultiplier)";
            case TestableType::Blocked:           return             "Cache-aware-blocked                    ";
            case TestableType::BlockedMatMul:     return             "Cache-aware-blocked  (MatrixMultiplier)";
            case TestableType::Recursive:         return std::format("Recursive until size {}               ", val); break;
            case TestableType::Strassen:          return std::format("Strassen  until size {}               ", val); break;
            case TestableType::Hybrid:            return             "Hybrid               (MatrixMultiplier)";
            case TestableType::Multithreaded:     return             "Multithreaded hybrid (MatrixMultiplier)";
        }
        return "";
    }

    std::string short_name, name;
    std::variant<
        MultiplierType,
        std::function<MultiplierType(int)>,
        std::function<MultiplierType(int, int, int)>
    > f;
};

struct TestConfig
{
    inline static std::set<TestableType> default_tests = {
        {TestableType::BetterNaive},
        {TestableType::BetterNaiveMatMul},
        {TestableType::Blocked},
        {TestableType::BlockedMatMul},
        {TestableType::Recursive, 32},
        {TestableType::Recursive, 64},
        {TestableType::Strassen, 32},
        {TestableType::Strassen, 64},
        {TestableType::Hybrid},
        {TestableType::Multithreaded}
    };
    inline static std::set<TestableType> all_tests = {
        {TestableType::Naive},
        {TestableType::NaiveMatMul},
        {TestableType::BetterNaive},
        {TestableType::BetterNaiveMatMul},
        {TestableType::Blocked},
        {TestableType::BlockedMatMul},
        {TestableType::Recursive, 16},
        {TestableType::Recursive, 32},
        {TestableType::Recursive, 64},
        {TestableType::Recursive, 128},
        {TestableType::Strassen, 16},
        {TestableType::Strassen, 32},
        {TestableType::Strassen, 64},
        {TestableType::Strassen, 128},
        {TestableType::Hybrid},
        {TestableType::Multithreaded}
    };

    inline static std::set<std::array<int, 3>> default_sizes = {
        {16,   16,   16},
        {32,   32,   32},
        {64,   64,   64},
        {128,  128,  128},
        {256,  256,  256},
        {512,  512,  512},
        {1024, 1024, 1024}    
    };

    std::vector<Testable> tests;
    std::set<std::array<int, 3>> sizes;
    bool verify_results = false;

    std::set<TestableType> tests_to_run;
    bool mult_parse_mode_with = true;

    bool size_parse_mode_with = true;

    void set_starting_mults(std::string_view arg)
    {
        if (arg == "default")
        {
            tests_to_run = default_tests;
            return;
        }
        if (arg == "all")
        {
            tests_to_run = all_tests;
            return;
        }
    }

    void set_starting_sizes_default()
    {
        sizes = default_sizes;
    }

    static std::optional<TestableType> parse_multiplier(std::string_view arg)
    {
        if (arg == "naive")               return {TestableType::Naive};
        if (arg == "naive_MatMul")        return {TestableType::NaiveMatMul};
        if (arg == "better_naive")        return {TestableType::BetterNaive};
        if (arg == "better_naive_MatMul") return {TestableType::BetterNaiveMatMul};
        if (arg == "blocked")             return {TestableType::Blocked};
        if (arg == "blocked_MatMul")      return {TestableType::BlockedMatMul};
        if (arg == "hybrid")              return {TestableType::Hybrid};
        if (arg == "multithreaded")       return {TestableType::Multithreaded};

        if (arg.starts_with("recursive") || arg.starts_with("strassen"))
        {
            std::string_view name = arg.starts_with("recursive") ? "recursive" : "strassen";
            int size;
            try
            {
                size = std::stoi(std::string(arg.substr(name.size())));
            }
            catch (std::invalid_argument& e)
            {
                std::cerr << "Invalid size format {" << arg << "} skipped\n";
                return std::nullopt;
            }
            if (size <= 0)
            {
                std::cerr << "Argument with negative size {" << arg << "} skipped\n";
                return std::nullopt;
            }
            if (name == "recursive") return TestableType{TestableType::Recursive, size};
            if (name == "strassen")  return TestableType{TestableType::Strassen,  size};
        }

        return std::nullopt;
    }

    static std::optional<std::array<int, 3>> parse_size(std::string_view arg)
    {
        auto sizes = arg | std::views::split('_');
        if (std::distance(sizes.begin(), sizes.end()) != 3)
        {
            std::cerr << "Invalid size format {" << arg << "} skipped\n";
            return std::nullopt;
        }
        std::array<int, 3> size{};
        try
        {
            std::ranges::copy(sizes | std::views::transform([](auto&& s){ return std::stoi(std::string(s.begin(), s.end())); }), size.begin());
        }
        catch(std::invalid_argument& e)
        {
            std::cerr << "Invalid size format {" << arg << "} skipped\n";
            return std::nullopt;
        }
        if (size[0] <= 0 || size[1] <= 0 || size[2] <= 0)
        {
            std::cerr << "Negative size {" << arg << "} skipped\n";
            return std::nullopt;
        }
        
        return size;
    } 

    void parse_config(zen::cmd_args& args)
    {
        verify_results = args.is_present("--verify");

        for (bool first = true; const auto& arg: args.get_options("--mult")) 
        {
            if (first && (arg == "default" || arg == "all")) 
            {
                first = false;
                set_starting_mults(arg);
                continue;
            } 
            if (arg == "with") 
            {
                mult_parse_mode_with = true;
                continue;
            }
            if (arg == "without") 
            {
                mult_parse_mode_with = false;
                continue;
            }
            if (auto type = parse_multiplier(arg)) 
            {
                if (mult_parse_mode_with) tests_to_run.insert(*type);
                else                      tests_to_run.erase(*type);
            }
        }
        for (const auto& arg: tests_to_run) tests.emplace_back(arg);

        for (bool first = true; const auto& arg: args.get_options("--sizes")) 
        {
            if (first && arg == "default") 
            {
                first = false;
                set_starting_sizes_default();
                continue;
            } 
            if (arg == "with") 
            {
                size_parse_mode_with = true;
                continue;
            }
            if (arg == "without") 
            {
                size_parse_mode_with = false;
                continue;
            }
            if (auto size = parse_size(arg)) 
            {
                if (size_parse_mode_with) sizes.insert(*size);
                else                      sizes.erase(*size);
            }
        }
    }
};


int main(int argc, char* argv[])
{
    zen::cmd_args args(argv, argc);

    if (argc == 1) 
    {
        std::cout << "Usage: " << argv[0] << " [--help | -h] [--verify] [--mult [default | all] [[with | without] <test_name>]...] [--sizes [default] [[with | without] <size1_size2_size3>]...]\n";
        std::cout << "Use --help or -h for detailed instructions.\n";
        return 0;
    }
    if (args.is_present("-h") || args.is_present("--help"))
    {
        std::cout << "Usage: " << args.first() << " [--verify] [--mult [default | all] [[with | without] <test_name>]...] [--sizes [default] [[with | without] <size1_size2_size3>]...]\n";
        std::cout << "Available multipliers:\n";
        for (const auto& type: 
            {
                TestableType::Naive,
                TestableType::NaiveMatMul,
                TestableType::BetterNaive,
                TestableType::BetterNaiveMatMul,
                TestableType::Blocked,
                TestableType::BlockedMatMul,
                TestableType::Recursive,
                TestableType::Strassen,
                TestableType::Hybrid,
                TestableType::Multithreaded
            })
        {
            auto name = Testable::short_name_from_type(type);
            if (name == "recursive" || name == "strassen")
                std::cout << "\t" << name << "<size>: ";
            else std::cout << "\t" << name << ": ";
            std::cout << Testable::name_from_type(type, std::optional{"<size>"}) << "\n";
        }
        std::cout << "\n\tdefault:\n";
        for (const auto& type: TestConfig::default_tests)
            std::cout << "\t\t" << Testable::short_name_from_type(type) << "\n";
        
        std::cout << "\n\tall:\n";
        for (const auto& type: TestConfig::all_tests)
            std::cout << "\t\t" << Testable::short_name_from_type(type) << "\n";

        std::cout << "\ndefault sizes:\n";
        for (const auto& size: TestConfig::default_sizes)
            std::cout << "\t" << size[0] << '_' << size[1] << '_' << size[2] << "\n";

        return 0;
    }

    TestConfig config;
    config.parse_config(args);

    std::cout << "CTEST_FULL_OUTPUT\n";
    
    for (const auto& size: config.sizes)
    {
        int N = size[0];
        int M = size[1];
        int P = size[2];
        std::cout << "N: " << N << ", M: " << M << ", P: " << P << '\n';
        for (const auto& test: config.tests)
        {
            if (std::holds_alternative<Testable::MultiplierType>(test.f))
            {
                auto f = std::get<Testable::MultiplierType>(test.f);
                print_test(test.name, f, N, M, P, config.verify_results);
            }
            else if (std::holds_alternative<std::function<Testable::MultiplierType(int, int, int)>>(test.f))
            {
                auto f = std::get<std::function<Testable::MultiplierType(int, int, int)>>(test.f)(N, M, P);
                print_test(test.name, f, N, M, P, config.verify_results);
            }
        }
    }

    return 0;
}