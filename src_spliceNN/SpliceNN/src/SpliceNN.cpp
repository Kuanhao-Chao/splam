#include <torch/torch.h>
#include <torch/script.h>
#include <iostream>
#include <memory>
#include <gclib/GArgs.h>

#define VERSION "0.0.6"

using namespace std;

const char* USAGE="SpliceNN v" VERSION "\n"
"==================\n"
"The TieCov utility can take the output file produced by TieBrush and generate the following auxiliary files:\n"
" 1. BedGraph file with the coverage data\n"
" 2. Junction BED file\n"
" 3. a heatmap BED that uses color intensity to represent the number of samples that contain each position\n"
"==================\n"
"\n"
" usage: tiecov [-s out.sample] [-c out.coverage] [-j out.junctions] [-W] input\n"
"\n"
" Input arguments (required): \n"
"  input\t\talignment file in SAM/BAM/CRAM format\n"
"       "
"\n"
" Optional arguments (at least one of -s/-c/-j must be specified):\n"
"  -h,--help\tShow this help message and exit\n"
"  --version\tShow program version and exit\n"
"  -s\t\tBedGraph file with an estimate of the number of samples\n"
"    \t\twhich contain alignments for each interval.\n"
"  -c\t\tBedGraph (or BedWig with '-W') file with coverage\n"
"    \t\tfor all mapped bases.\n"
"  -j\t\tBED file with coverage of all splice-junctions\n"
"    \t\tin the input file.\n"
"  -W\t\tsave coverage in BigWig format. Default output\n"
"    \t\tis in Bed format\n";


int main(int argc, const char* argv[]) {
  if (argc != 2) {
    cerr << "You must provide the path to the model!" << endl;
    return -1;    
  }

  torch::Tensor tensor = torch::rand({2, 3});
  std::cout << tensor << std::endl;
  const char *banner = R"""(
  ███████╗██████╗ ██╗      █████╗ ███╗   ███╗██╗
  ██╔════╝██╔══██╗██║     ██╔══██╗████╗ ████║██║
  ███████╗██████╔╝██║     ███████║██╔████╔██║██║
  ╚════██║██╔═══╝ ██║     ██╔══██║██║╚██╔╝██║╚═╝
  ███████║██║     ███████╗██║  ██║██║ ╚═╝ ██║██╗
  ╚══════╝╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝
  )""";
  std::cout << banner << std::endl;
  
  torch::jit::script::Module module;

  module = torch::jit::load(argv[1]);

  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    std::cout << "Loading "<< argv[1]<< std::endl;
    module = torch::jit::load(argv[1]);
  }
  catch (const c10::Error& e) {
    std::cerr << "error loading the model\n";
    return -1;
  }
  std::cout << "Model "<< argv[1]<<" loaded fine\n";
}

void processOptions(int argc, char* argv[]) {
    GArgs args(argc, argv, "help;verbose;version;DVWhc:s:j:");
    args.printError(USAGE, true);
    if (args.getOpt('h') || args.getOpt("help")) {
        GMessage(USAGE);
        exit(1);
    }

    if (args.getOpt("version")) {
        fprintf(stdout,"%s\n", VERSION);
        exit(0);
    }

    if ((args.getOpt('c') || args.getOpt('s') || args.getOpt('j'))==0){
        GMessage(USAGE);
        GMessage("\nError: at least one of -c/-j/-s arguments required!\n");
        exit(1);
    }

    verbose=(args.getOpt("verbose")!=NULL || args.getOpt('V')!=NULL);
    bigwig=args.getOpt('W')!=NULL;

    if (verbose) {
        fprintf(stderr, "Running TieCov " VERSION ". Command line:\n");
        args.printCmdLine(stderr);
    }
    covfname=args.getOpt('c');
    jfname=args.getOpt('j');
    sfname=args.getOpt('s');

    covfname_bw=args.getOpt('c');
    jfname_bw=args.getOpt('j');
    sfname_bw=args.getOpt('s');

    if (args.startNonOpt()==0) {
        GMessage(USAGE);
        GMessage("\nError: no input file provided!\n");
        exit(1);
    }
    
    infname=args.nextNonOpt();
}