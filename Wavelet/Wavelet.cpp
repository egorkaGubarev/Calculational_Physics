#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef unsigned int uint;
typedef float my_type;

template<typename type>
type f(const type t)
{
    const type pi = (type) std::acos(-1);
    const type result = (2 + 8 * std::cos(7 * 2 * pi * t) - 10 * std::sin(2 * pi * t)) / 10;
    return result;
}

template<typename type>
std::vector<type> get_upper_level_scaling_coeffs(const uint max_level, type (*f)(const type))
{
    const uint coeffs_amount = (uint) std::pow(2, max_level);
    std::vector<type> coeffs(coeffs_amount);

    for (uint coeff_number = 0; coeff_number < coeffs_amount; ++ coeff_number) {
        const type coeff = (*f)((type)coeff_number / coeffs_amount);
        coeffs[coeff_number] = coeff;
    }

    return coeffs;
}

template<typename type>
void get_scaling_coeffs(std::vector<std::vector<type>>& levels, const uint max_level)
{
    const type h0 = (type) 0.4829629131445341;
    const type h1 = (type) 0.8365163037378079;
    const type h2 = (type) 0.2241438680420133;
    const type h3 = (type) -0.1294095225512603;

    for (int level_number = max_level - 1; level_number >= 0; --level_number) {
        const uint upper_level_number = level_number + 1;
        const uint coeffs_amount = (uint) std::pow(2, level_number);

        const std::vector<type>& upper_coefs = levels[upper_level_number];
        std::vector<type>& coeffs = levels[level_number];
        
        coeffs.resize(coeffs_amount);

        for (uint coeff_number = 0; coeff_number < coeffs_amount; ++coeff_number) {
            const uint upper_coeffs_amount = (uint) std::pow(2, upper_level_number);
            const uint last_upper_coeff_index = upper_coeffs_amount - 1;
            const uint last_coeff_index = coeffs_amount - 1;

            type coeff_0 = 0, coeff_1 = 0, coeff_2 = 0, coeff_3 = 0;

            coeff_1 = upper_coefs[coeff_number * 2];
            coeff_2 = upper_coefs[coeff_number * 2 + 1];

            if (coeff_number == 0) {
                coeff_0 = upper_coefs[last_upper_coeff_index];
            }
            else {
                coeff_0 = upper_coefs[coeff_number * 2 - 1];
            }

            if (coeff_number == last_coeff_index) {
                coeff_3 = upper_coefs[0];
            }
            else {
                coeff_3 = upper_coefs[coeff_number * 2 + 2];
            }

            const type coeff = (type) (coeff_0 * h0 + coeff_1 * h1 + coeff_2 * h2 + coeff_3 * h3);
            coeffs[coeff_number] = coeff;
        }
    }
}

template<typename type>
void get_wavelet_coeffs(std::vector<std::vector<type>>& levels, std::vector<std::vector<my_type>>& wavelet_data, const uint max_level)
{
    const type h0 = (type) 0.4829629131445341;
    const type h1 = (type) 0.8365163037378079;
    const type h2 = (type) 0.2241438680420133;
    const type h3 = (type) -0.1294095225512603;

    const type g0 = h3;
    const type g1 = -h2;
    const type g2 = h1;
    const type g3 = -h0;

    for (int level_number = max_level - 1; level_number >= 0; --level_number) {
        const uint upper_level_number = level_number + 1;
        const uint coeffs_amount = (uint)std::pow(2, level_number);

        const std::vector<type>& upper_coefs = levels[upper_level_number];
        std::vector<type>& coeffs = wavelet_data[level_number];

        coeffs.resize(coeffs_amount);

        for (uint coeff_number = 0; coeff_number < coeffs_amount; ++coeff_number) {
            const uint upper_coeffs_amount = (uint)std::pow(2, upper_level_number);
            const uint last_upper_coeff_index = upper_coeffs_amount - 1;
            const uint last_coeff_index = coeffs_amount - 1;

            type coeff_0 = 0, coeff_1 = 0, coeff_2 = 0, coeff_3 = 0;

            coeff_1 = upper_coefs[coeff_number * 2];
            coeff_2 = upper_coefs[coeff_number * 2 + 1];

            if (coeff_number == 0) {
                coeff_0 = upper_coefs[last_upper_coeff_index];
            }
            else {
                coeff_0 = upper_coefs[coeff_number * 2 - 1];
            }

            if (coeff_number == last_coeff_index) {
                coeff_3 = upper_coefs[0];
            }
            else {
                coeff_3 = upper_coefs[coeff_number * 2 + 2];
            }

            const type coeff = (type)(coeff_0 * g0 + coeff_1 * g1 + coeff_2 * g2 + coeff_3 * g3);
            coeffs[coeff_number] = coeff;
        }
    }
}

template<typename type>
void save(const std::vector<std::vector<type>>& data, std::ofstream& file)
{
    for (const std::vector<type>& level : data) {

        for (const type value : level) {
            file << value << ' ';
        }

        file << '\n';
    }
}

int main()
{
    const uint max_level = 10;
    const uint levels_amount = max_level + 1;

    const uint resolution = (uint) std::pow(2, max_level);

    std::vector<std::vector<my_type>> scaling_data(levels_amount);
    std::vector<std::vector<my_type>> wavelet_data(max_level);

    std::string scaling_file_name = "scaling.txt";
    std::string wavelet_file_name = "wavelet.txt";

    std::ofstream scaling, wavelet;

    my_type(*f_ptr)(const my_type) = nullptr;
    f_ptr = &f;

    std::vector<my_type> data = get_upper_level_scaling_coeffs<my_type>(max_level, f_ptr);
    scaling_data[max_level] = data;

    get_scaling_coeffs<my_type>(scaling_data, max_level);
    get_wavelet_coeffs<my_type>(scaling_data, wavelet_data, max_level);

    scaling.open(scaling_file_name);
    wavelet.open(wavelet_file_name);

    save(scaling_data, scaling);
    save(wavelet_data, wavelet);

    scaling.close();
    wavelet.close();

    return 0;
}
