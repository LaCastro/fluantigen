/* Implements an individual-based model in which the infection's genealogical history is tracked through time */

public class Mutantigen {
    public static void main(String[] args) {

		// initialize random number generator
		cern.jet.random.AbstractDistribution.makeDefaultGenerator();
		
		String parametersFile = args[1];
		// initialize static parameters
		Parameters.load(parametersFile);		
		Parameters.initialize();
	    
                
		// initialize antigenic tree
        AntigenicTree.initialize();
               
		String dirFileName = args[0];
		String dirMainOutput = args[2];
		// run simulation
		Simulation sim = new Simulation();
		sim.run(dirFileName, dirMainOutput);	
		
		
	}
   	
}
