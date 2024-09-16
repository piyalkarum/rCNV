#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::CharacterVector process_snp(Rcpp::CharacterVector X, bool AD) {
  int n = X.size();
  std::vector<std::vector<double>> gg(n);

  // Split the string and convert to numeric
  for (int i = 0; i < n; ++i) {
    std::string cell = Rcpp::as<std::string>(X[i]);
    std::stringstream ss(cell);
    std::string item;
    while (std::getline(ss, item, ',')) {
      gg[i].push_back(std::stod(item));
    }
  }

  // Calculate column means
  std::vector<double> col_means(gg[0].size(), 0.0);
  for (const auto& row : gg) {
    for (size_t k = 0; k < row.size(); k++) {
      col_means[k] += row[k];
    }
  }
  for (auto& mean : col_means) {
    mean /= n;
  }

  // Remove the column with the minimum mean
  auto min_col_it = std::min_element(col_means.begin(), col_means.end());
  size_t min_col_idx = std::distance(col_means.begin(), min_col_it);

  for (auto& row : gg) {
    row.erase(row.begin() + min_col_idx);
  }

  // Convert back to required format
  Rcpp::CharacterVector result(n);
  for (int i = 0; i < n; ++i) {
    if (AD) {
      result[i] = std::to_string((int)gg[i][0]) + "," + std::to_string((int)gg[i][1]);
    } else {
      result[i] = std::to_string(gg[i][0]) + "," + std::to_string(gg[i][1]);
    }
  }

  return result;
}

// hetTgen function Helper function to split a string based on a delimiter
std::vector<std::string> split(const std::string &str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// [[Rcpp::export]]
CharacterVector process_column1(CharacterVector column, std::string info_type, int max_adn, IntegerVector ind_adn) {
    int n = column.size();
    CharacterVector result(n);

    for (int row = 0; row < n; ++row) {
        std::string cell = as<std::string>(column[row]);
        std::vector<std::string> tmp = split(cell, ':');

        if (tmp.size() <= static_cast<size_t>(ind_adn[row] - 1)) {
            result[row] = ""; // Handle case where ind_adn[row] is out of bounds
            continue;
        }

        if (info_type == "AD-tot") {
            std::vector<std::string> adn_parts = split(tmp[ind_adn[row] - 1], ',');
            if (adn_parts.size() == 2) {
                try {
                    double part1 = std::stod(adn_parts[0]);
                    double part2 = std::stod(adn_parts[1]);
                    result[row] = std::to_string(part1 + part2);
                } catch (const std::invalid_argument& e) {
                    result[row] = "0"; // Handle conversion error
                }
            } else {
                result[row] = "0"; // Handle unexpected format
            }
        } else {
            std::string adn_value = tmp[ind_adn[row] - 1];
            if (info_type != "DP" && (adn_value.empty() || adn_value == ".,.")) {
                adn_value = "./.";
            }
            result[row] = adn_value;
        }
    }

    return result;
}



