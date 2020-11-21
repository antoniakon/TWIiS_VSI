VariableSelectionMarch19

This project implements Bayesian variable selection, using the Kuo and Mallick approach, on the interactions of a hierarchical two-way Anova model.

There are four implementations:
- Both main and interaction effects are treated as asymmetric
- Both main and interaction effects are treated as symmetric
- Only the interaction effects are considered to be symmetric
- Only the main effects are considered to be symmetric

There are three branches: .
The first contains the fastest serial implementation. 
The branch "functionalSerial" is related to the serial implementation that includes more functional characteristics and is slower than the implementation in the "master".
The branch "parallel" contains the parallel version of the code based on the "functionalSerial" implementation that is concurrent safe.

To run the project after cloning the repository, open a terminal, change directory to the directory of the project and type:
git checkout NAME_OF_Branch
where NAME_OF_Branch is one of the: "master", "parallel" and "functionalSerial".

Then launch the sbt (Scala Build Tool) by typing:
java -jar PATH_TO_sbt-launch.jar

e.g. java -jar /usr/share/sbt/bin/sbt-launch.jar

To run the file with configuration type "run" and set the following parameters:
noOfIters: Int, //Number of iterations
thin: Int, //Thinning factor
burnIn: Int, //Burn-in period, implementationally it must be greater or equal to the thinning
logLikFlag: Boolean, //TRUE if we want to calculate and save the log-likelihood, e.g. to use it for the DIC criterion, FALSE otherwise
caseToRun: String, //One out of the four cases: "AsymmetricBoth", "SymmetricBoth", "SymmetricInteractions", "SymmetricMain"
pathToFiles: String, //Path to files 
inputFile: String, //csv input data file
outputFile: String, // csv file to save the output-sampling results
outputTimeFile: String //txt file to save the runtime

For example:
Iterations: 1,000,000, thinning: 100, burn-in period: 1,000

- For both main and interaction effects asymmetric
run 1000000 100 1000 true AsymmetricBoth ./data/15x20/ simulInterAsymmetricBoth.csv asymmetricBoth15x201m.csv asymmetricBoth15x201mRuntime.txt

- For both main and interaction effects symmetric
run 1000000 100 1000 true SymmetricBoth ./data/15x20/ simulInterSymmetricBoth.csv symmetricBoth15x201m.csv symmetricBoth15x201mRuntime.txt

- For only interaction effects as symmetric
run 1000000 100 1000 true SymmetricInteractions ./data/15x20/ simulInterSymmetricInters.csv symmetricInteractions15x201m.csv symmetricInteractions15x201mRuntime.txt

- For only main effects as symmetric
run 1000000 100 1000 true SymmetricMain ./data/15x20/ simulInterSymmetricMain.csv symmetricMain15x201m.csv symmetricMain15x201mRuntime.txt

Default parameters for model assumptions e.g. prior distributions are located in the file: "MainRunner.scala"