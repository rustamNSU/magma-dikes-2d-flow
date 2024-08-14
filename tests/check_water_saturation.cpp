#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>
#include "exprtk.hpp"

using SymbolTable = exprtk::symbol_table<double>;
using Expression  = exprtk::expression<double>;
using Parser      = exprtk::parser<double>;

const std::string expr_str = "sqrt(x) + 2x^2";

int main(int argc, char ** argv){
    double x = 0.0;
    SymbolTable symbol_table;
    symbol_table.add_variable("x", x);
    Expression expression;
    expression.register_symbol_table(symbol_table);
    Parser parser;
    bool is_compile = parser.compile(expr_str, expression);
    
    for (x = 0.0; x <= 1.0; x += 0.01){
        std::cout << x << " " << expression.value() << "\n";
    }
    return 0;
}