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
    const type result = (2 + 8 * std::cos(7 * 2 * pi * t) + 9 * std::sin(2 * pi * t));
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

            type coeff_0 = upper_coefs[(coeff_number * 2) % upper_coeffs_amount] ;
            type coeff_1 = upper_coefs[(coeff_number * 2 + 1) % upper_coeffs_amount];
            type coeff_2 = upper_coefs[(coeff_number * 2 + 2) % upper_coeffs_amount];
            type coeff_3 = upper_coefs[(coeff_number * 2 + 3) % upper_coeffs_amount];

            
            const type coeff = (type) (coeff_0 * h0 + coeff_1 * h1 + coeff_2 * h2 + coeff_3 * h3);
            coeffs[coeff_number] = coeff;
        }
    }
}

template<typename type>
void get_single_wavelet(std::vector<std::vector<type>>& scaling_coeffs)
{
    const type lowest_scaling = 0;
    const size_t levels = std::size(scaling_coeffs);
    const size_t wavelet_levels = levels - 1;
    std::vector<std::vector<type>> wavelet_coeffs(wavelet_levels);
    
    // Parsing wavelet coeffs table
    for (size_t i = 0; i < wavelet_levels; ++i) {
        const size_t coeffs_amount = (size_t) std::pow(2, i);
        std::vector<type>& level = wavelet_coeffs[i];
        level.resize(coeffs_amount);
    }

    std::vector<type>& lowest_wavelet = wavelet_coeffs[0];
    lowest_wavelet[0] = 1;

    make_inverse_transform<type>(lowest_scaling, wavelet_coeffs, scaling_coeffs);
}

template<typename type>
void get_wavelet_coeffs(std::vector<std::vector<type>>& levels, std::vector<std::vector<type>>& wavelet_data, const uint max_level)
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
            type coeff_0 = upper_coefs[(coeff_number * 2) % upper_coeffs_amount];

            type coeff_1 = upper_coefs[(coeff_number * 2 + 1) % upper_coeffs_amount];
            type coeff_2 = upper_coefs[(coeff_number * 2 + 2) % upper_coeffs_amount];
            type coeff_3 = upper_coefs[(coeff_number * 2 + 3) % upper_coeffs_amount];

            const type coeff = (type)(coeff_0 * g0 + coeff_1 * g1 + coeff_2 * g2 + coeff_3 * g3);
            coeffs[coeff_number] = coeff;
        }
    }
}

template<typename type>
void make_inverse_transform(const type lowest_scaling_coeff, const std::vector<std::vector<type>>& wavelet_coeffs,
    std::vector<std::vector<type>>& scaling_coeffs)
{
    const type h0 = (type)0.4829629131445341;
    const type h1 = (type)0.8365163037378079;
    const type h2 = (type)0.2241438680420133;
    const type h3 = (type)-0.1294095225512603;

    const type g0 = h3;
    const type g1 = -h2;
    const type g2 = h1;
    const type g3 = -h0;

    const size_t levels = std::size(wavelet_coeffs) + 1;

    scaling_coeffs.resize(levels);
    std::vector<type> lowest{lowest_scaling_coeff};
    scaling_coeffs[0] = lowest;

    // Itterating through all levels up
    for (size_t level_number = 1; level_number < levels; ++level_number) {
        const size_t prev_level_number = level_number - 1;

        const size_t coeffs_amount = (size_t) std::pow(2, level_number);
        size_t lower_coeffs_amount;

        // Correct handling of 0^2 = 1
        if (prev_level_number == 0) {
            lower_coeffs_amount = 1;
        }
        else {
            lower_coeffs_amount = (size_t)std::pow(2, prev_level_number);
        }
        
        std::vector<type> coeffs(coeffs_amount);

        // Getting all scaling coeffs of the level
        for (size_t coeff_number = 0; coeff_number < coeffs_amount; ++coeff_number) {
            const std::vector<type>& lower_scaling = scaling_coeffs[prev_level_number];
            const std::vector<type>& lower_wavelet = wavelet_coeffs[prev_level_number];

            // Chose of coeffs depends on odd or even number
            if (coeff_number % 2 == 0) {
                const int i = (coeff_number - 2) / 2 % lower_coeffs_amount;

                const type coeff_h_2 = lower_scaling[std::abs(i)];
                const type coeff_g_2 = lower_wavelet[std::abs(i)];
                const type coeff_h_0 = lower_scaling[coeff_number / 2];
                const type coeff_g_0 = lower_wavelet[coeff_number / 2];

                const type coeff = coeff_h_0 * h0 + coeff_h_2 * h2 + coeff_g_0 * g0 + coeff_g_2 * g2;
                coeffs[coeff_number] = coeff;
            }
            else {
                const int i = (coeff_number - 2) / 2 % lower_coeffs_amount;

                const type coeff_h_3 = lower_scaling[std::abs(i)];
                const type coeff_g_3 = lower_wavelet[std::abs(i)];
                const type coeff_h_1 = lower_scaling[coeff_number / 2];
                const type coeff_g_1 = lower_wavelet[coeff_number / 2];

                const type coeff = coeff_h_1 * h1 + coeff_h_3 * h3 + coeff_g_1 * g1 + coeff_g_3 * g3;
                coeffs[coeff_number] = coeff;
            }
        }
        scaling_coeffs[level_number] = coeffs;
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
    std::vector<std::vector<my_type>> scaling_data_of_single(levels_amount);
    std::vector<std::vector<my_type>> wavelet_data(max_level);
    std::vector<std::vector<my_type>> resored_scaling_coeffs;

    std::string scaling_file_name = "scaling.txt";
    std::string restored_file_name = "restore.txt";
    std::string wavelet_file_name = "wavelet.txt";
    std::string single_wavelet_file_name = "single.txt";

    std::ofstream scaling, wavelet, restore, single;

    my_type(*f_ptr)(const my_type) = nullptr;
    f_ptr = &f;

    std::vector<my_type> data = get_upper_level_scaling_coeffs<my_type>(max_level, f_ptr);
    scaling_data[max_level] = data;

    get_scaling_coeffs<my_type>(scaling_data, max_level);
    get_wavelet_coeffs<my_type>(scaling_data, wavelet_data, max_level);

    const std::vector<my_type>& lowest_scaling_level = scaling_data[0];
    const my_type lowest_scale_coeff = lowest_scaling_level[0];
    make_inverse_transform<my_type>(lowest_scale_coeff, wavelet_data, resored_scaling_coeffs);
    get_single_wavelet(scaling_data_of_single);

    scaling.open(scaling_file_name);
    wavelet.open(wavelet_file_name);
    restore.open(restored_file_name);
    single.open(single_wavelet_file_name);

    save(scaling_data, scaling);
    save(wavelet_data, wavelet);
    save(resored_scaling_coeffs, restore);
    save(scaling_data_of_single, single);

    scaling.close();
    wavelet.close();
    restore.close();
    single.close();

    return 0;
}
