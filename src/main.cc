/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "main.h"

// boost components
#include <boost/graph/iteration_macros.hpp>


// global debug variable
bool debug = false;

//! main program starts here
int main(int ac, char* av[]) {

    // initialize command line options
    boost::program_options::variables_map vm = init_options(ac, av);
    bool debug = vm["debug"].as<bool>();
    bool verbose = vm["verbose"].as<bool>();
    
    unsigned int number_of_designs = 10;
    if (vm.count("num")) {
        number_of_designs = vm["num"].as<unsigned int>();
    }
    
    std::string mode = "sample";
    if (vm.count("mode")) {
        mode = vm["mode"].as<std::string>();
    }
    
    // initialize mersenne twister with our seed
    unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
    if (vm.count("seed")) {
        seed = vm["seed"].as<unsigned long>();
    }
    if (verbose) {
        std::cerr << "Using this seed: " << seed << std::endl;
    }
    design::initialize_library(debug);

    // input handling ( we read from std::in per default and switch to a file if it is given in the --in option
    // std::in will provide a pseudo interface to enter structures directly!
    std::vector<std::string> structures;

    if (vm.count("in")) {
        std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
        if (infile->is_open()) {
            structures = read_input(infile);
            infile->close();
            std::cerr << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
        } else {
            std::cerr << "Unable to open file";
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "Input structures (one per line); @ to quit" << std::endl
                << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
        structures = read_input(&std::cin); // read infile into array
    }

    std::ostream* out = &std::cout;
    if (vm.count("out")) {
        out = new std::ofstream(vm["out"].as<std::string>(), std::ofstream::out);
    }
    // find sequence constraints in structures
    std::string constraints = "";
    std::regex seq ("[ACGTUWSMKRYBDHVN\\-]{1,}");
    std::regex str ("[\\(\\)\\.]{1,}");
    
    for (auto s = structures.begin(); s != structures.end();) {
        if (std::regex_match (*s, seq)) {
            constraints = *s;
            s = structures.erase(s);
        } else if (!std::regex_match (*s, str)) {
            std::cerr << "Unknown characters in line: " << *s << std::endl;
            exit(EXIT_FAILURE);
        } else {
            s++;
        }
    }
    
    if (debug) {
        std::cerr << "structures: " << std::endl 
                << structures << std::endl
                << "constraints: " << std::endl
                << constraints << std::endl;
    }
    
    design::DependencyGraph<std::mt19937>* dependency_graph;
    try {
        dependency_graph = new design::DependencyGraph<std::mt19937>(structures, constraints, seed);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (verbose) {
        std::cerr << "Size of solution space: " << dependency_graph->number_of_sequences() << std::endl;
    }
    // get an initial sequence
    dependency_graph->set_sequence(); // color the graph and get the sequence
    
    while (number_of_designs > 0) {
        try {
            if (mode == "sample")
                dependency_graph->set_sequence(); // color the graph and get the sequence
            else if (mode == "mutate-global")
                dependency_graph->mutate_global(); 
            else if (mode == "mutate-local")
                dependency_graph->mutate_local();
            else
                dependency_graph->set_sequence();
        } catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
        
        std::string sequence = dependency_graph->get_sequence();
        *out << sequence << std::endl;
        number_of_designs--;
    }
    
    delete dependency_graph;
    return EXIT_SUCCESS;
    
}

std::vector<std::string> read_input(std::istream * in) {
    // read input file
    std::string line;
    std::vector<std::string> structures;
    while (!in->eof()) {
        getline(*in, line);
        if (line == "@") {
            std::fclose(stdin);
        } else if (line.length() != 0) {
            structures.push_back(line);
        }
    }

    // exit if there is no input
    if (structures.empty()) {
        exit(EXIT_FAILURE);
    }

    if (debug) {
        std::cerr << "Read following structures:" << std::endl;
    }
    // check if structures have equal length
    unsigned int length = 0;
    for (auto elem : structures) {
        if (debug) {
            std::cerr << elem << std::endl;
        }
        if ((length != elem.length()) && (length != 0)) {
            std::cerr << "Structures have unequal length." << std::endl;
            exit(EXIT_FAILURE);
        }
        length = elem.length();
    }
    return structures;
}

boost::program_options::variables_map init_options(int ac, char* av[]) {
    // boost option parser
    // http://www.boost.org/doc/libs/1_53_0/doc/html/program_options/tutorial.html
    namespace po = boost::program_options;
    // Group of options that will be allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "print help message")
            ("verbose,v", po::bool_switch()->default_value(false)->zero_tokens(), "be verbose")
            ("debug,d", po::bool_switch()->default_value(false)->zero_tokens(), "be verbose for debugging")
            ;

    // Group of options that will be allowed on command line and in a config file
    po::options_description config("Program options");
    config.add_options()
            ("in,i", po::value<std::string>(), "input file which contains the structures [string]")
            ("out,o", po::value<std::string>(), "output file for writing the sequences [string]")
            ("mode,m", po::value<std::string>(), "mode for sequence generation\n\tsampling: stochastic sampling of all positions\n\
                                                \tmutate-global: Sample a initial sequence and then only mutate one connected component\n\
                                                \tmutate-local: Mutate only paths starting from an initial sequence  [string]")
            ("seed,s", po::value<unsigned long>(), "random number generator seed [unsigned long]")
            ("num,n", po::value<unsigned int>(), "number of designs (default: 10) [unsigned int]")
            ;

    po::positional_options_description p;
    p.add("in", 1).add("out", 2);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cmdline_options << "\n";
        exit(EXIT_FAILURE);
    }

    return vm;
}

