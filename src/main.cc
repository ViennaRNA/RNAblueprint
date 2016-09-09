/* This program reads RNA secondary structures in dot-bracket notation as well as
 * sequence constraints in IUPAC code and fairly samples RNA sequences compatible
 * to both inputs.
 *
 * @date 25.03.2013
 * @author Stefan Hammer <s.hammer@univie.ac.at>
 * @copyright GPLv3
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
    unsigned int number_of_designs = vm["num"].as<unsigned int>();
    std::string mode = vm["mode"].as<std::string>();
    
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
    std::string constraints;
    std::string start_seq;

    if (vm.count("in")) {
        std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
        if (infile->is_open()) {
            std::tie(structures, constraints, start_seq) = read_input(infile);
            infile->close();
            std::cerr << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
        } else {
            std::cerr << "Unable to open file";
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "Input structures (one per line), constraints and start sequence; @ to quit" << std::endl
                << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
        std::tie(structures, constraints, start_seq) = read_input(&std::cin); // read infile into array
    }

    std::ostream* out = &std::cout;
    if (vm.count("out")) {
        out = new std::ofstream(vm["out"].as<std::string>(), std::ofstream::out);
    }
    
    design::DependencyGraph<std::mt19937>* dependency_graph;
    try {
        dependency_graph = new design::DependencyGraph<std::mt19937>(structures, constraints, seed);
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    if (vm.count("graphml")) {
        std::ostream* graphml = new std::ofstream(vm["graphml"].as<std::string>(), std::ofstream::out);
        *graphml << dependency_graph->get_graphml();
        delete graphml;
    }

    if (verbose) {
        std::cerr << "Size of solution space: " << dependency_graph->number_of_sequences() << std::endl;
    }
    // get or set an initial sequence
    if (start_seq == "") {
        dependency_graph->sample(); // color the graph and get the sequence
    } else {
        dependency_graph->set_sequence(start_seq);
    }
    
    while (number_of_designs > 0) {
        try {
            if (mode == "sample")
                dependency_graph->sample(); // color the graph and get the sequence
            else if (mode == "sample-global")
                dependency_graph->sample_global(); 
            else if (mode == "sample-local")
                dependency_graph->sample_local();
            else if (mode == "global-neighbors") {
                dependency_graph->revert_sequence();
                dependency_graph->sample_global();
            } else if (mode == "local-neighbors") {
                dependency_graph->revert_sequence();
                dependency_graph->sample_local();
            } else
                dependency_graph->sample();
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

std::tuple<std::vector<std::string>, std::string, std::string > read_input(std::istream * in) {
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
    
        // find sequence constraints in structures
    std::string constraints = "";
    std::string start_seq = "";
    
    for (auto s = structures.begin(); s != structures.end();) {
        // find illegal character positions
        std::size_t illegal_constraint = s->find_first_not_of("ACGTUWSMKRYBDHVN-&+");
        std::size_t illegal_structure = s->find_first_not_of("().[]{}<>&+");
        // check if it is a constraint
        if (illegal_constraint == std::string::npos) {
            // check if it is a start sequence
            std::size_t illegal_sequence = s->find_first_not_of("ACGTU&+");
            if (illegal_sequence == std::string::npos && start_seq == "")
                start_seq = *s;
            else if (constraints == "")
                constraints = *s;
            else {
                std::cerr << "Too many constraints or start sequences: " << *s << std::endl;
                exit(EXIT_FAILURE);
            }
            s = structures.erase(s);
        // catch if there is a illegal character
        } else if (illegal_structure != std::string::npos) {
            std::cerr << "Illegal character [" << (*s)[illegal_structure] << "] in structure: " << *s << std::endl;
            exit(EXIT_FAILURE);
        } else {
            s++;
        }
    }
    
    if (debug) {
        std::cerr << "structures: " << std::endl 
                << structures << std::endl
                << "constraints: " << std::endl
                << constraints << std::endl
                << "start_seq: " << std::endl
                << start_seq << std::endl;
    }
    
    return std::make_tuple(structures, constraints, start_seq);
}

boost::program_options::variables_map init_options(int ac, char* av[]) {
    // boost option parser
    // http://www.boost.org/doc/libs/1_53_0/doc/html/program_options/tutorial.html
    namespace po = boost::program_options;
    po::options_description desc(     
        "This program reads RNA secondary structures in dot-bracket notation as well as\n"
        "sequence constraints in IUPAC code and fairly samples RNA sequences compatible\n"
        "to both inputs");
    // Group of options that will be allowed only on command line
    po::options_description generic("Generic Options");
    generic.add_options()
            ("help,h", "print help message")
            ("verbose,v", po::bool_switch()->default_value(false)->zero_tokens(), "be verbose")
            ("debug,d", po::bool_switch()->default_value(false)->zero_tokens(), "be verbose for debugging")
            ;

    // Group of options that will be allowed on command line and in a config file
    po::options_description config("Program options");
    config.add_options()
            ("in,i", po::value<std::string>(), "input file which contains the structures, sequence constraints and the start sequence [string]\n"
                                    "structures: \tsecondary structures in dot-bracket notation. one structure per input line\n"
                                    "sequence constraints: \tPermanent sequence constraints in IUPAC notation [ACGTUWSMKRYBDHVN] (optional)\n"
                                    "start sequence: \t A initial RNA sequence to start the sampling from [ACGU] (optional)")
            ("out,o", po::value<std::string>(), "output file for writing the sequences (default: stdout) [string]")
            ("graphml,g", po::value<std::string>(), "write a GraphML file representing the dependency graph to the given filename (optional) [string]")
            ("mode,m", po::value<std::string>()->default_value("sampling"), "mode for sequence generation [string]:\n"
                                    "sampling: \tstochastic sampling of all positions (default)\n"
                                    "sample-global: \tOnly sample one connected component at a time starting from an initial sequence\n"
                                    "sample-local: \tSample only single paths starting from an initial sequence\n"
                                    "global-neighbors: \tOnly find neighboring sequences to the initial start sequence by sampling one connected component only\n"
                                    "local-neighbors: \tOnly find neighboring sequences to the initial start sequence by sampling one path only\n"
            )
            ("seed,s", po::value<unsigned long>(), "random number generator seed [unsigned long]")
            ("num,n", po::value<unsigned int>()->default_value(10), "number of designs (default: 10) [unsigned int]")
            ;

    po::positional_options_description p;
    p.add("in", 1).add("out", 2);

    po::options_description cmdline_options;
    cmdline_options.add(desc).add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cmdline_options << "\n";
        exit(EXIT_FAILURE);
    }

    return vm;
}

