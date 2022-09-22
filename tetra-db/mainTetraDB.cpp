#include "tetra-db.hpp"

int main(){
    std::string starTable = "star_table.bin";
    std::string catalog = "catalog.bin";
    std::string properties = "properties.bin";
    std::cout << "hello" << std::endl;
    Tetra::GenerateDatabase(12, starTable, catalog, properties);
    return 0;
}