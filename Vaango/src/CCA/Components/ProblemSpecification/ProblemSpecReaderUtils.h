/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Limited, New
 * Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <Core/Exceptions/ProblemSetupException.h>

#include <algorithm>
#include <iomanip>
#include <string>

#include <libxml/parser.h>
#include <libxml/tree.h>

namespace Vaango {
namespace ProblemSpec {

// Helper for 'collapse' to remove all whitespace
void
collapse(std::string& s)
{
  s.erase(std::remove_if(
            s.begin(),
            s.end(),
            [](char c) { return std::isspace(static_cast<unsigned char>(c)); }),
          s.end());
}

// Helper for 'replace_substring'
void
replace_substring(std::string& s,
                  const std::string& from,
                  const std::string& to)
{
  size_t start_pos = 0;
  while ((start_pos = s.find(from, start_pos)) != std::string::npos) {
    s.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Handles case where 'to' contains 'from'
  }
}

std::string
get_error_info(const xmlNode* node)
{
  std::ostringstream error;

  if (node->_private == nullptr) {
    // All nodes of the ups_spec.xml will have _private set, but nodes coming
    // from the .ups file being validated may not.  However, they will have a
    // doc pointer that has the file information.
    //
    // Both ATTRIBUTE and TEXT node line numbers aren't part of those
    // type of nodes... we have to look at their parents to get the
    // real value.  (This occurs because we cast all nodes into
    // generic xmlNodes before they get to this portion of the code,
    // and this is why we need to do this test.)
    if (node->type == XML_ATTRIBUTE_NODE || node->type == XML_TEXT_NODE) {
      error << "See file: " << (const char*)(node->doc->URL) << " (line #"
            << node->parent->line << ")";
    } else {
      error << "See file: " << (const char*)(node->doc->URL) << " (line #"
            << node->line << ")";
    }
  } else {
    std::string file = *(std::string*)(node->_private);
    error << "See file: " << file << " (line #" << node->line << ")";
  }

  return error.str();
}

// valid_values may be:  "positive" | "*" | "num, num" which means min, max
// See ProblemSpecReader.h for more info on 'validValues'.
// An empty validValues means anything is valid.
void
validate_double(double value,
                const std::vector<std::string>& validValues,
                const std::string& completeName)
{
  if (validValues.empty()) { // Use .empty() for clarity
    return;
  }

  if (validValues.size() == 1) {
    // String comparison is fine here
    if (validValues[0] == "positive") {
      if (value < 0) {
        std::ostringstream error;
        error << std::setprecision(12); // Keep precision as in original
        error << "<" << completeName << ">: Specified value '" << value
              << "' is not 'positive' (as required).";
        throw Uintah::ProblemSetupException(error.str(), __FILE__, __LINE__);
      }
    } else {
      // If there's only one validValue, but it's not "positive",
      // the original code implicitly falls through to the 'else' branch,
      // which throws "Invalid 'validValues' string.". We'll maintain that.
      throw Uintah::ProblemSetupException(
        completeName + " - Invalid 'validValues' string. Expected 'positive' "
                       "for single value validation.",
        __FILE__,
        __LINE__);
    }
  } else if (validValues.size() == 2) {
    double min, max;
    try {
      size_t pos;
      min = std::stod(validValues[0], &pos); // Convert first string to min
      // Verify no trailing non-whitespace characters
      for (size_t i = pos; i < validValues[0].length(); ++i) {
        if (!std::isspace(static_cast<unsigned char>(validValues[0][i]))) {
          throw std::invalid_argument("trailing characters");
        }
      }

      max = std::stod(validValues[1], &pos); // Convert second string to max
      // Verify no trailing non-whitespace characters
      for (size_t i = pos; i < validValues[1].length(); ++i) {
        if (!std::isspace(static_cast<unsigned char>(validValues[1][i]))) {
          throw std::invalid_argument("trailing characters");
        }
      }

    } catch (const std::invalid_argument& e) {
      // Catches if conversion fails (e.g., "abc") or if trailing chars exist
      std::ostringstream error;
      error << "<" << completeName << ">: Invalid min/max range values. "
            << "Expected two double values, but found non-numeric data or "
               "trailing characters in '"
            << validValues[0] << "' or '" << validValues[1] << "'. "
            << "Details: " << e.what();
      throw Uintah::ProblemSetupException(error.str(), __FILE__, __LINE__);
    } catch (const std::out_of_range& e) {
      // Catches if the number is too large/small for a double
      std::ostringstream error;
      error << "<" << completeName
            << ">: Min/max range values are out of double range. "
            << "Values: '" << validValues[0] << "', '" << validValues[1]
            << "'. "
            << "Details: " << e.what();
      throw Uintah::ProblemSetupException(error.str(), __FILE__, __LINE__);
    }

    if (value < min || value > max) {
      std::ostringstream error;
      error << std::setprecision(12); // Keep precision as in original
      error << "<" << completeName << "> - "
            << "Specified value '" << value << "' is outside of valid range ("
            << min << ", " << max << ").";
      throw Uintah::ProblemSetupException(error.str(), __FILE__, __LINE__);
    }
  } else {
    // This 'else' branch handles validValues.size() > 2
    throw Uintah::ProblemSetupException(
      completeName +
        " - Invalid 'validValues' string. Expected 0, 1, or 2 elements.",
      __FILE__,
      __LINE__);
  }
}

// Read string, convert, and validate a single double
void
read_and_validate_single_double(const std::string& text,
                                const xmlNode* node,
                                const std::string& classType,
                                const std::vector<std::string>& validValues,
                                const std::string& completeName)
{
  try {
    // std::stod converts string to double. It throws on invalid input.
    // It also ignores leading/trailing whitespace.
    size_t pos;
    double value = std::stod(text, &pos);

    // Check if there's any non-whitespace garbage after the number
    // This is equivalent to sscanf's "only a double exists in the text" check
    bool has_trailing_chars = false;
    for (size_t i = pos; i < text.length(); ++i) {
      if (!std::isspace(static_cast<unsigned char>(text[i]))) {
        has_trailing_chars = true;
        break;
      }
    }

    if (has_trailing_chars) {
      throw Uintah::ProblemSetupException(
        classType + " <" + completeName +
          "> should have only a double value (but has: '" + text +
          "'). Please fix XML in .ups file or correct validation Tag "
          "list.\n" +
          get_error_info(node),
        __FILE__,
        __LINE__);
    }

    validate_double(
      value, validValues, completeName); // Call existing custom validation

  } catch (const std::invalid_argument& e) {

    // Catches if no conversion could be performed (e.g., "abc")
    throw Uintah::ProblemSetupException(
      classType + " <" + completeName +
        "> should have a double value (but has: '" + text +
        "'). Please fix XML in .ups file or correct validation Tag "
        "list.\n" +
        get_error_info(node),
      __FILE__,
      __LINE__);

  } catch (const std::out_of_range& e) {

    // Catches if the converted value is out of the range of a double
    throw Uintah::ProblemSetupException(
      classType + " <" + completeName +
        "> double value is out of range (but has: '" + text +
        "'). Please fix XML in .ups file or correct validation Tag "
        "list.\n" +
        get_error_info(node),
      __FILE__,
      __LINE__);
  }
}

// Read string, convert, and validate a single integers
void
read_and_validate_single_integer(const std::string& text,
                                 const xmlNode* node,
                                 const std::string& classType,
                                 const std::vector<std::string>& validValues,
                                 const std::string& completeName)
{
  try {
    // std::stoi converts string to int. Throws on invalid input/out of range.
    // It also ignores leading/trailing whitespace.
    size_t pos;
    int value = std::stoi(text, &pos);

    // Check for trailing characters, similar to DOUBLE case
    bool has_trailing_chars = false;
    for (size_t i = pos; i < text.length(); ++i) {
      if (!std::isspace(static_cast<unsigned char>(text[i]))) {
        has_trailing_chars = true;
        break;
      }
    }

    if (has_trailing_chars) {
      throw Uintah::ProblemSetupException(
        classType + " <" + completeName +
          "> should have only an integer value (but has: '" + text +
          "'). Please fix XML in .ups file or correct validation Tag "
          "list.\n" +
          get_error_info(node),
        __FILE__,
        __LINE__);
    }

    validate_double(static_cast<double>(value),
                    validValues,
                    completeName); // Still calls validateDouble

  } catch (const std::invalid_argument& e) {
    throw Uintah::ProblemSetupException(
      classType + " <" + completeName +
        "> should have an integer value (but has: '" + text +
        "'). Please fix XML in .ups file or correct validation Tag "
        "list.\n" +
        get_error_info(node),
      __FILE__,
      __LINE__);
  } catch (const std::out_of_range& e) {
    throw Uintah::ProblemSetupException(
      classType + " <" + completeName +
        "> integer value is out of range (but has: '" + text +
        "'). Please fix XML in .ups file or correct validation Tag "
        "list.\n" +
        get_error_info(node),
      __FILE__,
      __LINE__);
  }
}

// Read string, convert, and validate multiple doubles
void
read_and_validate_multiple_doubles(const std::string& text,
                                   const xmlNode* node,
                                   const std::string& classType,
                                   const std::string& completeName)
{
  // Modern C++ splitting using std::string::find and std::string::substr
  // Delimiters: "[,]"
  std::string delimiters = "[,]";
  size_t lastPos         = 0;
  size_t pos             = text.find_first_of(delimiters, 0);

  while (true) {
    std::string token_str;
    if (pos == std::string::npos) {
      token_str = text.substr(lastPos);
    } else {
      token_str = text.substr(lastPos, pos - lastPos);
    }

    if (!token_str.empty()) {
      try {
        size_t token_pos;
        [[maybe_unused]] double result = std::stod(token_str, &token_pos);

        // Ensure no trailing characters in the individual token
        bool has_trailing_chars = false;
        for (size_t i = token_pos; i < token_str.length(); ++i) {
          if (!std::isspace(static_cast<unsigned char>(token_str[i]))) {
            has_trailing_chars = true;
            break;
          }
        }

        if (has_trailing_chars) {
          throw Uintah::ProblemSetupException(
            classType + " ('" + completeName +
              "') contains a non-double token ('" + token_str + "') in '" +
              text +
              "'. Please fix XML in .ups file or correct validation Tag "
              "list.\n" +
              get_error_info(node),
            __FILE__,
            __LINE__);
        }

        // No further validation for individual doubles in the original code
      } catch (const std::invalid_argument& e) {
        throw Uintah::ProblemSetupException(
          classType + " ('" + completeName +
            "') contains a non-double token ('" + token_str + "') in '" + text +
            "'. Please fix XML in .ups file or correct validation Tag "
            "list.\n" +
            get_error_info(node),
          __FILE__,
          __LINE__);
      } catch (const std::out_of_range& e) {
        throw Uintah::ProblemSetupException(
          classType + " ('" + completeName +
            "') contains an out-of-range double token ('" + token_str +
            "') in '" + text +
            "'. Please fix XML in .ups file or correct validation Tag "
            "list.\n" +
            get_error_info(node),
          __FILE__,
          __LINE__);
      }
    }

    if (pos == std::string::npos) {
      break;
    }
    lastPos = pos + 1;
    pos     = text.find_first_of(delimiters, lastPos);
  }
}

// Read string, convert, and validate multiple integers
void
read_and_validate_multiple_integers(const std::string& text,
                                    const xmlNode* node,
                                    const std::string& classType,
                                    const std::string& completeName)
{
  // Check for decimals as in original code
  if (text.find(".") != std::string::npos) { // Use std::string::npos
    throw Uintah::ProblemSetupException(
      classType + " ('" + completeName +
        "') should have multiple integer values (but has: '" + text +
        "'). Please fix XML in .ups file or correct validation Tag "
        "list.\n" +
        get_error_info(node),
      __FILE__,
      __LINE__);
  }

  // Modern C++ splitting using std::string::find and std::string::substr
  // Delimiters: "[,]"
  std::string delimiters = "[,]";
  size_t lastPos         = 0;
  size_t pos             = text.find_first_of(delimiters, 0);

  while (true) {
    std::string token_str;
    if (pos == std::string::npos) {
      token_str = text.substr(lastPos);
    } else {
      token_str = text.substr(lastPos, pos - lastPos);
    }

    // Skip empty tokens that might result from adjacent delimiters
    // or leading/trailing delimiters.
    if (!token_str.empty()) {
      try {
        size_t token_pos;
        [[maybe_unused]] int result = std::stoi(token_str, &token_pos);

        // Ensure no trailing characters in the individual token
        bool has_trailing_chars = false;
        for (size_t i = token_pos; i < token_str.length(); ++i) {
          if (!std::isspace(static_cast<unsigned char>(token_str[i]))) {
            has_trailing_chars = true;
            break;
          }
        }

        if (has_trailing_chars) {
          throw Uintah::ProblemSetupException(
            classType + " ('" + completeName +
              "') contains a non-integer token ('" + token_str + "') in '" +
              text +
              "'. Please fix XML in .ups file or correct validation Tag "
              "list.\n" +
              get_error_info(node),
            __FILE__,
            __LINE__);
        }

        // No further validation for individual integers in the original code,
        // so we just ensure successful conversion.
      } catch (const std::invalid_argument& e) {
        throw Uintah::ProblemSetupException(
          classType + " ('" + completeName +
            "') contains a non-integer token ('" + token_str + "') in '" +
            text +
            "'. Please fix XML in .ups file or correct validation Tag "
            "list.\n" +
            get_error_info(node),
          __FILE__,
          __LINE__);
      } catch (const std::out_of_range& e) {
        throw Uintah::ProblemSetupException(
          classType + " ('" + completeName +
            "') contains an out-of-range integer token ('" + token_str +
            "') in '" + text +
            "'. Please fix XML in .ups file or correct validation Tag "
            "list.\n" +
            get_error_info(node),
          __FILE__,
          __LINE__);
      }
    }

    if (pos == std::string::npos) {
      break; // No more delimiters, finished
    }
    lastPos = pos + 1;
    pos     = text.find_first_of(delimiters, lastPos);
  }
}

bool
validate_vector(const std::string& text)
{
    // Create a mutable copy of the input string
    std::string cleanText = text;

    // 1. Remove all whitespace characters from the string.
    // This is more robust than just replacing " " with ""
    cleanText.erase(std::remove_if(cleanText.begin(), cleanText.end(),
                                   [](unsigned char c){ return std::isspace(c); }),
                    cleanText.end());

    // 2. Initial format check: Minimum length, starts with '[', ends with ']'
    if (cleanText.length() < 3 || cleanText.front() != '[' || cleanText.back() != ']') {
        return false;
    }

    // 3. Count commas: A vector [x,y,z] needs exactly two commas.
    // std::count is efficient for single characters.
    int numCommas = std::count(cleanText.begin(), cleanText.end(), ',');
    if (numCommas != 2) {
        return false;
    }

    // 4. Use std::istringstream for structured parsing.
    // This allows reading characters and numbers in sequence.
    std::istringstream iss(cleanText);

    double val1, val2, val3;
    char open_bracket, comma1, comma2, close_bracket;

    // Attempt to extract the expected format: '[double,double,double]'
    iss >> open_bracket >> val1 >> comma1 >> val2 >> comma2 >> val3 >> close_bracket;

    // 5. Check if the extraction was successful and characters matched expectations.
    // iss.fail() indicates a conversion error (e.g., trying to read a double from "abc")
    // or a format mismatch (e.g., trying to read a char 'X' when it's 'Y').
    if (iss.fail() || iss.bad() ||
        open_bracket != '[' || comma1 != ',' || comma2 != ',' || close_bracket != ']')
    {
        return false;
    }

    // 6. Crucial: Check for any remaining characters in the stream after parsing the vector.
    // If !iss.eof() after consuming potential trailing whitespace, there's unparsed garbage.
    iss >> std::ws; // Consume any trailing whitespace
    if (!iss.eof()) {
        return false; // Trailing garbage found (e.g., "[1,2,3]ABC")
    }

    // If all checks pass, the string conforms to the expected vector format.
    return true;
}

// Read string, convert, and validate multiple integers
void
read_and_validate_multiple_vectors(const std::string& text,
                                   const xmlNode* node,
                                   const std::string& classType,
                                   const std::string& completeName)
{
  // Keep original logic for replace_substring, collapse, and substr checks
  // as they are already using std::string methods.
  std::string tempText = text;
  replace_substring(tempText, " ", "");
  replace_substring(tempText, "\t", "");
  collapse(tempText); // Remove all whitespace

  // Vector of Vectors starts with [[ and ends with ]]... verify this...
  if (tempText.length() < 4 || tempText.substr(0, 2) != "[[" ||
      tempText.substr(tempText.length() - 2, 2) != "]]") {
    throw Uintah::ProblemSetupException(
      "This does not look like a Vector of Vectors. "
      "Expected to find [[ and ]] but have this : '" + tempText +
        "'",
      __FILE__,
      __LINE__);
  }

  // Start after '[[', first vector starts after pos 1 in tempText.
  // original started at 1, but then substr was (pos, len)
  // so if text is "[[1,2,3],[4,5,6]]", first '[' is at index 2
  // so start at 2 to find ']' and then '['.
  unsigned int pos = 2; 

  // The original logic `pos = tempText.find("[", pos + 1);` for finding the
  // next vector's opening bracket is still good. The `while (pos <
  // tempText.length() - 2)` check is also fine.
  while (pos < tempText.length() - 2 &&
         pos !=
           std::string::npos) { // Added pos != std::string::npos for safety

    size_t end_bracket_pos = tempText.find("]", pos);
    if (end_bracket_pos == std::string::npos) {
      // Malformed string, no closing bracket found for a vector
      throw Uintah::ProblemSetupException(
        classType + " ('" + completeName +
          "') has a malformed MULTIPLE_VECTOR value (missing closing ']'). " +
          "Problematic part: '" + tempText.substr(pos) + "'. " +
          get_error_info(node),
        __FILE__,
        __LINE__);
    }

    std::string vectorStr = tempText.substr(pos, end_bracket_pos - pos + 1);

    if (!validate_vector(vectorStr)) { // Call existing custom validation
      throw Uintah::ProblemSetupException(
        classType + " ('" + completeName +
          "') should have a MULTIPLE_VECTOR value (but has: '" + text +
          "'). Please fix XML in .ups file or correct validation Tag "
          "list.\n" +
          get_error_info(node),
        __FILE__,
        __LINE__);
    }
    pos = tempText.find("[", pos + 1); // Find the next opening bracket
  }
}


} // end namespace ProblemSpec
} // end namespace Vaango