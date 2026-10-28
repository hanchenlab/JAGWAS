#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator> 
#include <numeric>
#include <string> 
#include <armadillo>
#include <cstdlib> 
#include <map>

// Dependencies for the log10P option
#include <pgamma.c>
#include <lgamma.c>
#include <fmax2.c>
#include <mlutils.c>
//#include <dpois.c>
#include <dnorm.c>
#include <pnorm.c>
#include <ftrunc.c>
#include <gamma.c>
#include <lgammacor.c>
#include <stirlerr.c>
#include <bd0.c>
#include <chebyshev.c>


// z-score input has the columns CHR SNP POS A1 A2 N AF1 BETA SE Zscore, so the z-score is the 10th column
void checkZscoreColumns(const std::vector<std::string>& data, const std::string& fileName) {
    if (data.size() < 10) {
        std::cerr << "Error: a row of " << fileName << " has " << data.size() << " columns. "
                  << "z-score input needs 10 columns (CHR SNP POS A1 A2 N AF1 BETA SE Zscore), with the z-score in column 10." << std::endl;
        std::exit(1);
    }
}


void processFiles(const std::string& outputPath, const std::string& Cor_M, int nrow , double MAF, bool score_test, bool Beta_se , bool logP, std::string delimiter, const std::vector<std::string>& fileNames) {
 
    int nline = nrow;
    double AF;
    double AF1;
    arma::mat cor;
    char delim;
    std::string value;
    std::string Value;
    std::string line;
    std::string Line;
    std::string data;
    std::size_t dim = fileNames.size();
    std::vector<std::ifstream> file_handles;
    std::ofstream results(outputPath);
    std::ostringstream oss;  //buffer
    if (logP == true) {
    oss << "CHR" << "\t" << "SNP" << "\t" << "POS" << "\t" << "A1" << "\t" << "A2" << "\t" << "N" << "\t" << "AF1" << "\t" << "log10_P"<< "\n";
    }
    else
    {
    oss << "CHR" << "\t" << "SNP" << "\t" << "POS" << "\t" << "A1" << "\t" << "A2" << "\t" << "N" << "\t" << "AF1" << "\t" << "P"<< "\n";
    }
    
    
//Set up delimiter
if ((delimiter[0] == '\\' && delimiter[1] == 't') || delimiter[0] == 't') {
        delim = '\t';
        }
        else if ((delimiter[0] == '\\' && delimiter[1] == '0') || delimiter[0] == '0') {
        delim = ' ';
        }
        else {
        delim = delimiter[0];
        }




    for (size_t f = 0; f < fileNames.size(); f++) {
        // std::ifstream currentFile = openFile(fileNames[f]);
        file_handles.emplace_back(fileNames[f]);
    }
   
    if (!cor.load(Cor_M, arma::auto_detect)) {
        std::cerr << "Error loading the matrix from file." << std::endl;
    }
    arma::mat cor_i = inv(cor);

    int p = 0;

//To skip the header
std::getline(file_handles[0], line);
for (size_t f = 1; f < fileNames.size(); f++) {
        std::getline(file_handles[f], Line);
    }



if (score_test == true) {

while (std::getline(file_handles[0], line)) {
  
         
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, delim)) data.push_back(value);
        AF = std::stod(data[4]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        double score = std::stod(data[5]);
        double var = std::stod(data[6]);
        zscores[0] = (score/var)/sqrt(1/var);
        
        bool skip_this_variant = false; // Flag to control skipping
// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Value;
        std::string Line;
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, delim)) data.push_back(Value);

        AF1 = std::stod(data[4]);
        if (AF1 != AF) {
                    std::cerr << "Warning: Skipping inconsistent variant at line " << p + 1
                              << ". Mismatch in file: " << fileNames[f] << std::endl;
                    skip_this_variant = true;

                    // IMPORTANT: Consume lines from remaining files to keep them in sync
                    for (size_t f_skip = f + 1; f_skip < fileNames.size(); ++f_skip) {
                        std::string dummyLine;
                        std::getline(file_handles[f_skip], dummyLine);
                    }
                    break; // Exit the inner for-loop
                }
        double score = std::stod(data[5]);
        double var = std::stod(data[6]);
        zscores[f] = (score/var)/sqrt(1/var);
    } 


    if (skip_this_variant) {
                p++;
                continue; // Skips to the next iteration of the main `while` loop
            }
    // --- This code only runs if the variant was fully consistent ---
    // write the variant columns only now, so a skipped variant leaves no partial line in the output
    oss << "NA" << "\t" << data[0] << "\t" << "NA" << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t";
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor_i * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2.0, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every n line 
            oss.str(std::string());
            oss.clear();
        }	
    }

    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
    }
    // do the matrix production, get the chi-sq test statistics and output the results  
   

        }    
    results << oss.str();    //for the remainder
    oss.str(std::string());
    oss.clear();
    results.close();

}







else {

    if (Beta_se == false) {
    while (std::getline(file_handles[0], line)) {
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, delim)) data.push_back(value);
        checkZscoreColumns(data, fileNames[0]);
        AF = std::stod(data[6]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        zscores[0] = std::stod(data[9]);
        
        bool skip_this_variant = false; // Flag to control skipping

// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, delim)) data.push_back(Value);
        checkZscoreColumns(data, fileNames[f]);
        AF1 = std::stod(data[6]);
        if (AF1 != AF) {
                    std::cerr << "Warning: Skipping inconsistent variant at line " << p + 1
                              << ". Mismatch in file: " << fileNames[f] << std::endl;
                    skip_this_variant = true;

                    // IMPORTANT: Consume lines from remaining files to keep them in sync
                    for (size_t f_skip = f + 1; f_skip < fileNames.size(); ++f_skip) {
                        std::string dummyLine;
                        std::getline(file_handles[f_skip], dummyLine);
                    }
                    break; // Exit the inner for-loop
                }
        zscores[f] = std::stod(data[9]);
    }
  


    if (skip_this_variant) {
                p++;
                continue; // Skips to the next iteration of the main `while` loop
            }
    // --- This code only runs if the variant was fully consistent ---
    // write the variant columns only now, so a skipped variant leaves no partial line in the output
    oss << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5] << "\t" << data[6] << "\t";
    // do the matrix production, get the chi-sq test statistics and output the results
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor_i * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2.0, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every nline 
            oss.str(std::string());
            oss.clear();
        }	
    }
    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
         }
    }
    }
       
    




if (Beta_se == true) {
    while (std::getline(file_handles[0], line)) {
  
         
        std::istringstream iss(line);
        std::vector<std::string> data;
        while (std::getline(iss, value, delim)) data.push_back(value);
        AF = std::stod(data[6]);
    if (AF>MAF && AF<(1-MAF))
    {
        std::vector<double> zscores(dim);
        double beta = std::stod(data[7]);
        double se = std::stod(data[8]);
        zscores[0] = beta/se;

        bool skip_this_variant = false; // Flag to control skipping

// chunk into different fields, when reading, save fields to string,  

       for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Value;
        std::string Line;
        std::getline(file_handles[f], Line);
        std::istringstream iss(Line);
        std::vector<std::string> data;
        while (std::getline(iss, Value, delim)) data.push_back(Value);

        AF1 = std::stod(data[6]);
        if (AF1 != AF) {
                    std::cerr << "Warning: Skipping inconsistent variant at line " << p + 1
                              << ". Mismatch in file: " << fileNames[f] << std::endl;
                    skip_this_variant = true;

                    // IMPORTANT: Consume lines from remaining files to keep them in sync
                    for (size_t f_skip = f + 1; f_skip < fileNames.size(); ++f_skip) {
                        std::string dummyLine;
                        std::getline(file_handles[f_skip], dummyLine);
                    }
                    break; // Exit the inner for-loop
                }
        double beta = std::stod(data[7]);
        double se = std::stod(data[8]);
        zscores[f] = beta/se; 
    } 

    if (skip_this_variant) {
                p++;
                continue; // Skips to the next iteration of the main `while` loop
            }
      
// --- This code only runs if the variant was fully consistent ---
    // write the variant columns only now, so a skipped variant leaves no partial line in the output
    oss << data[0] << "\t" << data[1] << "\t" << data[2] << "\t" << data[3] << "\t" << data[4] << "\t" << data[5] << "\t" << data[6] << "\t";
   arma::vec Zscores = zscores;
   arma::mat z = arma::reshape(Zscores, 1, Zscores.size());
   //arma::mat zt = arma::reshape(Zscores, Zscores.size(), 1);
   arma::mat zt = trans(z);
   arma::mat chisq = z * cor_i * zt;
   double chisqs = chisq(0,0);
   double pValue = pgamma(chisqs/4, dim/2.0, 0.5, 0, logP);
   if (logP == true) {
   pValue = pValue*(1/log(10));      
   }
   oss << pValue << "\n";

   p++;
    
        if (p % nline == 0) {
            results << oss.str(); // dump oss to results, every n line 
            oss.str(std::string());
            oss.clear();
        }	
    }

    else {
        for (size_t f = 1; f < fileNames.size(); f++) {
        std::string Line;
        std::getline(file_handles[f], Line);
    }
    }
    // do the matrix production, get the chi-sq test statistics and output the results  
   

        }    

    }

    results << oss.str();    //for the remainder
    oss.str(std::string());
    oss.clear();
    results.close();
}
}

        

int main(int argc, char* argv[]) {

    std::map<std::string, std::string> cmdArgs;
    std::map<std::string, bool> knownArgs {
        {"--outputFilePath", false},
        {"--cor_matrix", false},
        {"--nrow", false},
        {"--MAF", false},
        {"--score_test", false},
        {"--beta_se", false},
        {"--logP", false},
        {"--delim", false},
        {"--fileNames", false}
    };
    std::vector<std::string> fileNames;
    



for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg.substr(0, 2) == "--") {
            // Check if the argument is known
            if (knownArgs.find(arg) != knownArgs.end()) {
                if (arg == "--fileNames") {
                    // every value up to the next option is a file name
                    while (i + 1 < argc && std::string(argv[i + 1]).substr(0, 2) != "--") {
                        fileNames.push_back(argv[++i]);
                    }
                    if (fileNames.empty()) {
                        std::cerr << "Value for " << arg << " not specified." << std::endl;
                        return 1;
                    }
                    knownArgs[arg] = true;
                } else if (i + 1 < argc && std::string(argv[i + 1]).substr(0, 2) != "--") { // Make sure we are not at the end of argv!
                    cmdArgs[arg] = argv[i + 1];  // Add the argument to our map
                    knownArgs[arg] = true;  // Mark this argument as found
                    i++;  // Increment 'i' so we don't try to process the value as a key
                } else {
                    std::cerr << "Value for " << arg << " not specified." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Unknown argument: " << arg << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Unexpected value: " << arg << std::endl;
            return 1;
        }
    }

    bool missing = false;
    for (const auto& a : knownArgs) {
        if (!a.second) {
            std::cerr << "Missing argument: " << a.first << std::endl;
            missing = true;
        }
    }
    if (missing) {
        return 1;
    }




    // read the options by name, so they can be given in any order
    std::string outputFilePath(cmdArgs["--outputFilePath"]);
    std::string cor_m(cmdArgs["--cor_matrix"]);
    int nrow = std::atoi(cmdArgs["--nrow"].c_str());
    double MAF = std::stod(cmdArgs["--MAF"]);
    bool score_test = std::stoi(cmdArgs["--score_test"]) != 0;
    bool beta_se = std::stoi(cmdArgs["--beta_se"]) != 0;
    bool logP = std::stoi(cmdArgs["--logP"]) != 0;
    std::string delimiter(cmdArgs["--delim"]);
    

    processFiles(outputFilePath, cor_m, nrow,  MAF, score_test, beta_se, logP, delimiter, fileNames);
    return 0;

}



