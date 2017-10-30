/* Simulation functions, holds the host population */

import java.util.*;
import java.io.*;
import java.text.SimpleDateFormat;

//import com.javamex.classmexer.*;
import com.javamex.classmexer.*;

public class Simulation {

	// fields
	private String time;
	private List<HostPopulation> demes = new ArrayList<HostPopulation>();
	private double diversity;
	private double tmrca;
	private double netau;
	private double serialInterval;
	private double antigenicDiversity;
	private double meanLoad;

	private List<Double> diversityList = new ArrayList<Double>();
	private List<Double> tmrcaList = new ArrayList<Double>();	
	private List<Double> netauList = new ArrayList<Double>();		
	private List<Double> serialIntervalList = new ArrayList<Double>();		
	private List<Double> antigenicDiversityList = new ArrayList<Double>();
	private List<Double> nList = new ArrayList<Double>();
	private List<Double> sList = new ArrayList<Double>();	
	private List<Double> iList = new ArrayList<Double>();	
	private List<Double> rList = new ArrayList<Double>();		
	private List<Double> casesList = new ArrayList<Double>();

	private ArrayList<Double> daysList = new ArrayList<Double>();

	// Added these to track viral fitness dists
	private double meanRDist;
	private double varRDist;
	private double meanBetaDist;
	private double varBetaDist;
	private double meanSigmaDist;
	private double varSigmaDist;
	private double covBetaSigmaDist;


	// constructor
	public Simulation() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			if (Parameters.restartFromCheckpoint) {
				HostPopulation hp = new HostPopulation(i, true);
				demes.add(hp);
			}
			else {
				HostPopulation hp = new HostPopulation(i);
				demes.add(hp);
			}
		}
	}

	// methods

	public int getN() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getN();
		}
		return count;
	}

	public int getS() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getS();
		}
		return count;
	}	

	public int getI() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getI();
		}
		return count;
	}	

	public int getR() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getR();
		}
		return count;
	}		

	public int getCases() {
		int count = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			count += hp.getCases();
		}
		return count;
	}	

	public double getDiversity() {
		return diversity;
	}		

	public double getNetau() {
		return netau;
	}	

	public double getTmrca() {
		return tmrca;
	}	

	public double getSerialInterval() {
		return serialInterval;	
	}		

	public double getAntigenicDiversity() {
		return antigenicDiversity;
	}		

	// proportional to infecteds in each deme
	public int getRandomDeme() {
		int n = Random.nextInt(0,getN()-1);
		int d = 0;
		int target = (demes.get(0)).getN();
		while (n < target) {
			d += 1;
			target += (demes.get(d)).getN();
		}
		return d;
	}

	// return random virus proportional to worldwide prevalence
	public Virus getRandomInfection() {

		Virus v = null;

		if (getI() > 0) {

			// get deme proportional to prevalence
			int n = Random.nextInt(0,getI()-1);
			int d = 0;
			int target = (demes.get(0)).getI();
			while (d < Parameters.demeCount) {
				if (n < target) {
					break;
				} else {
					d++;
					target += (demes.get(d)).getI();
				}	
			}
			HostPopulation hp = demes.get(d);

			// return random infection from this deme
			if (hp.getI()>0) {
				Host h = hp.getRandomHostI();
				v = h.getInfection();
			}

		}

		return v;

	}

	public double getMeanLoad() {
		HostPopulation hp = demes.get(0);
		double mean = hp.getMeanLoad();
		return mean;
	}

	public int getAntigenicTypesCount() {
		HostPopulation hp = demes.get(0);
		int types = hp.getAntigenicCount();
		return types;
	}

	public ArrayList<Integer> getAntigenicTypes() {
		ArrayList<Integer> typeList = new ArrayList<Integer>();
		HostPopulation hp = demes.get(0);

		if (getI() > 0) {
			for (int i = 0; i < getI(); i++) {
				Host h = hp.getInfecteds().get(i);
				Virus v = h.getInfection();
				Phenotype p = v.getPhenotype();
				int type = p.antigenicType();
				if (!typeList.contains(type)) {
					typeList.add(type);
				}
			}
		}
		return typeList;
	}



	public int getCumlAntigenicTypes() {
		//int cumlTypes = AntigenicTree.nodeCount();
		ArrayList<Integer> treeTypes = AntigenicTree.getTreeTypes();
		int cumlTypes = treeTypes.size();
		return cumlTypes;
	}

	//        public void reconstructAntDynamics(ArrayList<Integer> types) {
	//            
	//            //Print to file
	//        }


	// return random host from random deme
	public Host getRandomHost() {
		int d = Random.nextInt(0,Parameters.demeCount-1);
		HostPopulation hp = demes.get(d);
		return hp.getRandomHost();
	}

	public double getAverageRisk(Phenotype p) {

		double averageRisk = 0;
		for (int i = 0; i < 10000; i++) {
			Host h = getRandomHost();
			Phenotype[] history = h.getHistory();
			averageRisk += p.riskOfInfection(history);
		}
		averageRisk /= 10000.0;
		return averageRisk;

	}

	public void printImmunity() {

		try {
			File immunityFile = new File("out.immunity");
			immunityFile.delete();
			immunityFile.createNewFile();
			PrintStream immunityStream = new PrintStream(immunityFile);

			for (double x = VirusTree.xMin; x <= VirusTree.xMax; x += 0.5) {
				for (double y = VirusTree.yMin; y <= VirusTree.yMax; y += 0.5) {

					Phenotype p = PhenotypeFactory.makeArbitaryPhenotype(x,y);
					double risk = getAverageRisk(p);
					immunityStream.printf("%.4f,", risk);

				}
				immunityStream.println();
			}

			immunityStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}

	}

	public void printHostPopulation() {

		try {
			File hostFile = new File("out.hosts");
			hostFile.delete();
			hostFile.createNewFile();
			PrintStream hostStream = new PrintStream(hostFile);
			for (int i = 0; i < Parameters.demeCount; i++) {
				HostPopulation hp = demes.get(i);
				hp.printHostPopulation(hostStream);
			}
			hostStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}

	}	

	public void makeTrunk() {
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.makeTrunk();
		}
	}	

	public void printState() {

		System.out.printf("%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\n", 
				Parameters.day, getDiversity(), getTmrca(),  getNetau(), getSerialInterval(), getAntigenicDiversity(), 
				getN(), getS(), getI(), getR(), getCases(), getMeanLoad(), getAntigenicTypesCount(), getCumlAntigenicTypes());

		if (Parameters.memoryProfiling && Parameters.day % 100 == 0) {
			long noBytes = 0;
			//long noBytes = MemoryUtil.deepMemoryUsageOf(this);
			//System.out.println("Total: " + noBytes);
			//HostPopulation hp = demes.get(1);
			//noBytes = MemoryUtil.deepMemoryUsageOf(hp);
			//System.out.println("One host population: " + noBytes);
			//HostPopulation hp = demes.get(0);
			//Host h = hp.getRandomHostS();
			//long noBytes = MemoryUtil.deepMemoryUsageOf(h);
			//System.out.println("One susceptible host with " +  h.getHistoryLength() + " previous infection: " + noBytes);
			//h.printHistory();
			if (getI() > 0) {
				//Virus v = getRandomInfection();
				//noBytes = MemoryUtil.memoryUsageOf(v);
				//System.out.println("One virus: " + noBytes);
				noBytes = MemoryUtil.deepMemoryUsageOf(VirusTree.getTips());
				double doubleBytes = (double) noBytes;
				double megaBytes = 6000.0;
				double conversion = 1050000.0;
				double percent = doubleBytes / (megaBytes * conversion);
				System.out.println("Virus tree: " + noBytes + " " + "(" + percent + ")");

				noBytes = MemoryUtil.deepMemoryUsageOf(AntigenicTree.getDistArray());
				doubleBytes = (double) noBytes;
				percent = doubleBytes / (megaBytes * conversion);
				System.out.println("Antigenic tree: " + noBytes + " " + "(" + percent + ")");
			}
		}

	}

	public void printHeader(PrintStream stream) {
		stream.print("date\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\ttotalN\ttotalS\ttotalI\ttotalR\ttotalCases");
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printHeader(stream);
		}
		stream.println();
	}

	public void printState(PrintStream stream) {
		stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", Parameters.getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printState(stream);
		}
		stream.println();
	}

	public void printMutations(PrintStream stream) {
		//stream.printf("%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d", Parameters.getDate(), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), getN(), getS(), getI(), getR(), getCases());
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.printMutations(stream);
		}
		stream.println();
	}

	public void printFitness(PrintStream stream) {
		stream.printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", Parameters.day-1, Parameters.getDate(), meanRDist, varRDist, meanBetaDist, varBetaDist, meanSigmaDist, varSigmaDist, covBetaSigmaDist);
		stream.println();
	}

	
	
	public void printFrequencies(PrintStream stream, ArrayList<Double> currentFrequencies, ArrayList<Integer> typeList) {
		// want to put current infected count here 
		int infected = this.getI();
		for (int i = 0; i < currentFrequencies.size(); i++) {
			stream.printf("%d\t%d\t%.6f\t%d", (Parameters.day-1), typeList.get(i), currentFrequencies.get(i), infected ); 
			stream.println();
		}
	}

	
	public void printTypeMutations(PrintStream stream, ArrayList<Double> typeMutations, ArrayList<Double> varMutationLoads, ArrayList<Integer> typeList){
		for (int i = 0; i < typeMutations.size(); i++) {
			stream.printf("%d\t%d\t%.6f\t%.6f", (Parameters.day-1), typeList.get(i), typeMutations.get(i), varMutationLoads.get(i)); 
			stream.println();
		}
	}
	
	
	public void printViralFitnessDistTypes(PrintStream stream, VirusFitnessDistType fitDistType) {
		stream.printf("%d\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f", (Parameters.day-1), Parameters.getDate(),
				fitDistType.antigenicType, fitDistType.meanMut, fitDistType.varMut, fitDistType.meanR, fitDistType.varR,
				fitDistType.meanBeta, fitDistType.varBeta, fitDistType.meanSigma, fitDistType.varSigma);
		stream.println();
	}
	
	public void printTypeImmunities(PrintStream stream, ArrayList<Double> typeImmunities, ArrayList<Integer> typeList) {
		for (int i = 0; i <typeImmunities.size(); i++) {
			stream.printf("%d\t%d\t%.6f", (Parameters.day-1), typeList.get(i), typeImmunities.get(i));
			stream.println();
		}
	}
	
	public void printTrackAntigens(PrintStream stream) {
		stream.printf("%d\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d", 
				(Parameters.day-1), getDiversity(), getTmrca(), getNetau(), getSerialInterval(), getAntigenicDiversity(), 
				getN(), getS(), getI(), getR(), getCases(), getMeanLoad(), getAntigenicTypesCount(),getCumlAntigenicTypes());
		stream.println();
	}

	public void printFitnessSamples(String time) {

		try {
			String fileName = time + "/out.fitnessSamples_t" + Double.toString(Parameters.day);
			File samplesFile = new File(fileName);
			samplesFile.delete();
			samplesFile.createNewFile();
			PrintStream samplesStream = new PrintStream(samplesFile);

			samplesStream.print("betas\tS/N\tR");

			// Get samples
			HostPopulation hp = demes.get(0);
			ArrayList<ArrayList<Double>> samplesData = hp.getViralFitnessSamples();

			int samples = samplesData.get(0).size();
			for (int s = 0; s < samples; s++) {
				double beta = samplesData.get(0).get(s);
				double sigma = samplesData.get(1).get(s);
				double R = samplesData.get(2).get(s);
				samplesStream.printf("%.6f\t%.6f\t%.6f", beta, sigma, R);
				samplesStream.println();
			}

			samplesStream.close();

		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}

	}

	public void printFitnessHeader(PrintStream stream) {
		stream.print("day\tsimDay\tmeanR\tvarR\tmeanBeta\tvarBeta\tmeanSigma\tvarSigma\tcovBetaSigma");
		stream.println();
	}
	
	public void printFitnessTypeHeader(PrintStream stream) {
		stream.print("day\tsimDay\tantigenType\tmeanMut\tvarMut\tmeanR\tvarR\tmeanBeta\tvarBeta\tmeanSigma\tvarSigma");
		stream.println();
	}

	public void printTrackAntigenHeader(PrintStream stream) {
		stream.print("day\tdiversity\ttmrca\tnetau\tserialInterval\tantigenicDiversity\tN\tS\tI\tR\tcases\tmeanLoad\tantigenicTypes\tcumulativeTypes");
		stream.println();
	}

	public void printTrackFrequenciesHeader(PrintStream stream) {
		stream.print("day\tantigentype\tfrequency\tinfected");
		stream.println();
	}
	
	public void printTypeMutationsHeader(PrintStream stream) {
		stream.print("day\tantigentype\tnumMutations\tvarMutations");
		stream.println();
	}

	public void printTypeImmunitiesHeader(PrintStream stream) {
		stream.print("day\tantigentype\timmunity");
		stream.println();;
	}
	
	public void printAntigenicDistances(ArrayList<Integer> types) {
		AntigenicTree.printDistanceMatrix(types, time);
	}

	public void printAntigenicDistancesBetweenTips() {
		AntigenicTree.printDistanceMatrixBetweenTips(time);
	}

	public void printSummary(String time) {

		try {
			String summaryName = time.concat("/out.summary.txt");
			File summaryFile = new File(summaryName);
			summaryFile.delete();
			summaryFile.createNewFile();
			PrintStream summaryStream = new PrintStream(summaryFile);
			summaryStream.printf("parameter\tfull\n");

			summaryStream.printf("diversity\t%.4f\n", mean(diversityList));
			summaryStream.printf("tmrca\t%.4f\n", mean(tmrcaList));
			summaryStream.printf("netau\t%.4f\n", mean(netauList));		
			summaryStream.printf("serialInterval\t%.5f\n", mean(serialIntervalList));	
			summaryStream.printf("antigenicDiversity\t%.4f\n", mean(antigenicDiversityList));	
			summaryStream.printf("N\t%.4f\n", mean(nList));		
			summaryStream.printf("S\t%.4f\n", mean(sList));		
			summaryStream.printf("I\t%.4f\n", mean(iList));		
			summaryStream.printf("R\t%.4f\n", mean(rList));		
			summaryStream.printf("cases\t%.4f\n", mean(casesList));					

			summaryStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}

	}

	private double mean(List<Double> list) {
		double mean = 0;
		if(!list.isEmpty()) {
			for (Double item : list) {
				mean += (double) item;
			}
			mean /= (double) list.size();
		}
		return mean;
	}	

	public void updateDiversity() {

		diversity = 0.0;
		tmrca = 0.0;
		antigenicDiversity = 0.0;		
		netau = 0.0;
		serialInterval = 0.0;

		double coalCount = 0.0;	
		double coalOpp = 0.0;
		double coalWindow = Parameters.netauWindow / 365.0;
		int sampleCount = Parameters.diversitySamplingCount;

		for (int i = 0; i < sampleCount; i++) {
			Virus vA = getRandomInfection();
			Virus vB = getRandomInfection();
			if (vA != null && vB != null) {
				double dist = vA.distance(vB);
				diversity += dist;
				if (dist > tmrca) {
					tmrca = dist; // if the distance is greater, keep putting back the most recent common ancestor?
				}
				antigenicDiversity += vA.antigenicDistance(vB);
				coalOpp += coalWindow;
				coalCount += vA.coalescence(vB, coalWindow); // TODO look this up: is there a coalescence event within x amount of time?
				serialInterval += vA.serialInterval(); // birth of virus A to its parents birth
			}
		}	

		diversity /= (double) sampleCount;
		tmrca /= 2.0;	
		netau = coalOpp / coalCount;
		serialInterval /= (double) sampleCount;
		antigenicDiversity /= (double) sampleCount;			

	}

	public void updateFitnessDists() {

		HostPopulation hp = demes.get(0);
		//double mean = hp.getAntigenicLoadDistribution(); // is this used???
		VirusFitnessDist fitDist = hp.getViralFitnessDistribution(); // VirusFitnessDist holds variables
		meanRDist = fitDist.meanR;
		varRDist = fitDist.varR;
		meanBetaDist = fitDist.meanBeta;
		varBetaDist = fitDist.varBeta;
		meanSigmaDist = fitDist.meanSigma;
		varSigmaDist = fitDist.varSigma;
		covBetaSigmaDist = fitDist.covBetaSigma;

	}

	// LC 
	
	public void updateFitnessDistsTypes(PrintStream stream) {
		HostPopulation hp = demes.get(0);
		ArrayList<Integer> typeList = hp.getAntigenicTypes();
		ArrayList<ArrayList<Integer>> mutationLoadDistributions = hp.getVarMutationType(typeList);
		//ArrayList<ArrayList<Double>> sigmaLoadDistributions = hp.getDistTypeImmunities(typeList);
		ArrayList<Double> meanImmunities = hp.getMeanTypeImmunities(typeList);
		ArrayList<Double> varImmunities = hp.getVarTypeImmunities(typeList);
		
		for (int t = 0; t < typeList.size(); t ++) {
		
			ArrayList<Integer> mutationList = mutationLoadDistributions.get(t);
			double meanSigma = meanImmunities.get(t);
			double varSigma = varImmunities.get(t);
			//ArrayList<Double> sigmaList = sigmaLoadDistributions.get(t);
			int type = typeList.get(t);
			
		VirusFitnessDistType fitnessDistsTypes = hp.getViralFitnessDistrubtionTypes(meanSigma, varSigma, mutationList, type);
		printViralFitnessDistTypes(stream, fitnessDistsTypes);
		
		}
	}
	
	
	public ArrayList<Double> updateFrequencies() {
		HostPopulation hp = demes.get(0);
		ArrayList<Integer> typeList = hp.getAntigenicTypes();
		//ArrayList<Integer> typeCounts = hp.getAntigenicTypeCounts(typeList);
		// need to pull out the index of typeList and typeCounts 
		ArrayList<Double> typeFrequencies = hp.getAntigenicTypeCounts(typeList);


		for (int i = 0; i < typeList.size(); i++) {
			//System.out.println(typeList.get(i) + " " + typeFrequencies.get(i));
		}
		return typeFrequencies;
	}
	// LC 
	public ArrayList<Double> updateMeanTypeImmunities() {
		HostPopulation hp = demes.get(0);
		ArrayList<Integer> typeList = hp.getAntigenicTypes();
		ArrayList<Double> meanTypeImmunities = hp.getMeanTypeImmunities(typeList);
		return(meanTypeImmunities);
	}
	
	// LC 
	public ArrayList<Double> updateMutationLoads() {
		HostPopulation hp = demes.get(0);
		ArrayList<Integer> typeList = hp.getAntigenicTypes();
		ArrayList<Double> typeMutationLoad = hp.getMutationTypeCounts(typeList);
		return typeMutationLoad;
	}
	
	// LC 
//	public ArrayList<Double> updateVarMutationTypes() {
//		HostPopulation hp = demes.get(0);
//		ArrayList<Integer> typeList = hp.getAntigenicTypes();
//		ArrayList<Double> varMutationLoadTypes = hp.getVarMutationType(typeList);
//		return(varMutationLoadTypes);
				
//	}

	public void pushLists() {
		diversityList.add(diversity);
		tmrcaList.add(tmrca);
		netauList.add(netau);
		serialIntervalList.add(serialInterval);
		antigenicDiversityList.add(antigenicDiversity);
		nList.add((double) getN());
		sList.add((double) getS());
		iList.add((double) getI());
		rList.add((double) getR());
		casesList.add((double) getCases());		
	}
	
	// LC 

	public void resetCases() {
		for (int i = 0; i < Parameters.demeCount; i++) {	
			HostPopulation hp = demes.get(i);
			hp.resetCases();
		}
	}

	public void stepForward() { // stepping forward the simualation

		for (int i = 0; i < Parameters.demeCount; i++) {		
			HostPopulation hp = demes.get(i); //the host population is an attribute of the simulation

			//hp.stepForward(); // Trevor's method

			hp.stepForwardMutAtTrans(); // My method with mutation at transmission

			for (int j = 0; j < Parameters.demeCount; j++) {
				if (i != j) {
					HostPopulation hpOther = demes.get(j);
					hp.betweenDemeContact(hpOther);
				}
			}
		}

		Parameters.day++;

	}

	public void run(String dirFileName, String dirMainOutput) {

		try {

			
			//Date now = new Date();
			//SimpleDateFormat dateFormat = new SimpleDateFormat("MM-dd-yyyy_hh-mm-ss");
			//String timeFormat = dateFormat.format(now);
			time = dirMainOutput + "/" + dirFileName;
			//time = dirFileName;
			new File(time).mkdirs();
			//dir.mkdir();

			System.setOut(new PrintStream(new FileOutputStream(time.concat("/out.console.txt"))));

			String timeseriesName = time.concat("/out.timeseries.txt");

			File seriesFile = new File(timeseriesName);	// changed to a text file	

			//File seriesFile = new File("out.timeseries.txt");	// changed to a text file	
			seriesFile.delete();
			seriesFile.createNewFile();
			PrintStream seriesStream = new PrintStream(seriesFile);
			System.out.println("day\tsimDay\toriAntigen\tdistance\tmutLoad\tpostAntigen"); // console output for antigenic mutations
			printHeader(seriesStream);

			String mutationSeriesName = time.concat("/out.mutationSeries.txt");
			File mutationFile = new File(mutationSeriesName);// changed to a text file
			mutationFile.delete();
			mutationFile.createNewFile();
			PrintStream mutationStream = new PrintStream(mutationFile);

			String viralFitnessName = time.concat("/out.viralFitnessSeries.txt");
			File fitnessFile = new File(viralFitnessName); // changed to a text file
			fitnessFile.delete();
			fitnessFile.createNewFile();
			PrintStream fitnessStream = new PrintStream(fitnessFile);
			printFitnessHeader(fitnessStream);

			// Lauren adding to keep track of when antigenic mutations occur 

			String antigenSeriesName = time.concat("/out.trackAntigenSeries.txt");
			File trackAntigenFile = new File(antigenSeriesName);
			trackAntigenFile.delete();
			trackAntigenFile.createNewFile();
			PrintStream trackAntigenStream = new PrintStream(trackAntigenFile);
			printTrackAntigenHeader(trackAntigenStream);

		 
			String antigenFrequenciesName = time.concat("/out.antigenFrequencies.txt");
			File trackFrequenciesFile = new File(antigenFrequenciesName);
			trackFrequenciesFile.delete();
			trackFrequenciesFile.createNewFile();
			PrintStream trackFrequenciesStream = new PrintStream(trackFrequenciesFile);
			printTrackFrequenciesHeader(trackFrequenciesStream);

/*			String typeMutationsName = time.concat("/out.typeMutations.txt");
			File trackTypeMutationsFile = new File(typeMutationsName);
			trackTypeMutationsFile.delete();
			trackTypeMutationsFile.createNewFile();
			PrintStream typeMutationsStream = new PrintStream(trackTypeMutationsFile);
			printTypeMutationsHeader(typeMutationsStream);
	*/		
			
			String typeViralFitness = time.concat("/out.typeViralFitness.txt");
			File trackViralFitnessFile = new File(typeViralFitness);
			trackViralFitnessFile.delete();
			trackViralFitnessFile.createNewFile();
			PrintStream typeViralFitnessStream = new PrintStream(trackViralFitnessFile);
			printFitnessTypeHeader(typeViralFitnessStream);
			
			for (int i = 0; i < Parameters.endDay; i++) {

				stepForward();


				/*if(Parameters.novelAntigen == true && Parameters.day > Parameters.burnin) {
					//System.out.println("Novel antigen-time to update");
					//updateDiversity(); // make sure this isn't pushing anything 
					//updateFitnessDists(); // make sure this isn't pushing anything 
					printTrackAntigens(trackAntigenStream);
				
					ArrayList<Double> currentFrequencies = updateFrequencies(); // New method to update frequencies of each antigen
					ArrayList<Integer> typeList = getAntigenicTypes();
					
					printFrequencies(trackFrequenciesStream, currentFrequencies, typeList);
					printFitness(fitnessStream);
					
					updateFitnessDistsTypes(typeViralFitnessStream);
		/*			ArrayList<Double> currentMutationLoads = updateMutationLoads(); // New method to update average mutation load of each antigen
					ArrayList<Double> currentVarMutationLoads = updateVarMutationTypes();
					printTypeMutations(typeMutationsStream, currentMutationLoads, currentVarMutationLoads, typeList);				
					ArrayList<Double> currentTypeImmunities = updateMeanTypeImmunities();  // New method to update mean type immunities
					printTypeImmunities(typeImmunitiesStream, currentTypeImmunities, typeList);
	
				}
				/* Here I want to be able to UpdateDiversity, UpdateFitnessDists
				 * Print out printTrackAntigens(trackAntigenStream) I think this is taken care of 
				 
				}
				*/
				if (Parameters.day % Parameters.printStep == 0 ) {
					daysList.add(Parameters.getDate());
					updateDiversity();
					//printState();
					if (Parameters.day > Parameters.burnin) {
						printState(seriesStream);
						printMutations(mutationStream);	
						printTrackAntigens(trackAntigenStream);
						//	if(Parameters.novelAntigen == false) {
						updateFitnessDists(); // new for tracking viral fitness distributions
						printFitness(fitnessStream);
						ArrayList<Double> currentFrequencies = updateFrequencies(); // New method to update frequencies of each antigen
						ArrayList<Integer> typeList = getAntigenicTypes();
						printFrequencies(trackFrequenciesStream, currentFrequencies, typeList);
						updateFitnessDistsTypes(typeViralFitnessStream);
						//}
						/*				ArrayList<Double> currentMutationLoads = updateMutationLoads(); // New method to print out average mutation load of each antigen
						ArrayList<Double> currentVarMutationLoads = updateVarMutationTypes();
						printTypeMutations(typeMutationsStream, currentMutationLoads, currentVarMutationLoads, typeList);
						ArrayList<Double> currentTypeImmunities = updateMeanTypeImmunities();  // New method to update mean type immunities
						printTypeImmunities(typeImmunitiesStream, currentTypeImmunities, typeList);
						 */					
						
					}
					//resetCases(); change so it's annual attack 
					pushLists();
				}
				
				Parameters.novelAntigen = false;
				
				if (Parameters.day % 365 == 0 ) {
					resetCases();
				}

				if (Parameters.day % Parameters.printFitSamplesStep == 0) {
					// Print list of randomly sampled virus's beta, S/N and R values
					printFitnessSamples(time);
				}

				if (getI()==0) {
					if (Parameters.repeatSim) {
						reset();
						i = 0; 
						seriesFile.delete();
						seriesFile.createNewFile();
						seriesStream = new PrintStream(seriesFile);
						printHeader(seriesStream);
					} else {
						break;
					}
				}
			}

			seriesStream.close();
			mutationStream.close();
			fitnessStream.close();
		} catch(IOException ex) {
			System.out.println("Could not write to file"); 
			System.exit(0);
		}	


		// Summary
		printSummary(time);

		// tree reduction
		/*VirusTree.pruneTips();
		VirusTree.markTips();		

		// tree prep
		makeTrunk();
		VirusTree.fillBackward();                                       // assign parents children			
		VirusTree.sortChildrenByDescendants();                          // sort children by number of descendents
		VirusTree.setLayoutByDescendants();
		VirusTree.getMRCASeries(daysList, time);                        // prints out time of MRCA in tree over time
		VirusTree.streamline();                                         // collapses nodes with single descentdent in tree			


		// reconstruct dynamics of antigenic types in tree
		VirusTree.getTreeTypes();                                       // get all antigenic types in tree (including internals)
		VirusTree.reconstructAntDynamics(daysList, time);               // reconstruct the population dynamics of all antigenic types in tree
		printAntigenicDistances(VirusTree.treeTypes);                   // print antigenic distance matrix for antigenic types in tree

		// rotation
		if (Parameters.pcaSamples) {
			VirusTree.rotate();
			VirusTree.flip();
		}

		// tip and tree output
		VirusTree.printTips(time);			
		VirusTree.printBranches(time);	
		VirusTree.printNewick(time);                                        // Print basic newick tree
		VirusTree.printSimMapLoad(time);                                    // Print SimMap tree with annotation for mutation load along lineages
		//VirusTree.printSimMapClades();
		VirusTree.printSimMapAntigenic(time);                               // Print SimMap tree with annotation for antigenic type along lineages
		VirusTree.printTrunkAntigenicShifts(time);                          // Print size of all antigenic transitions along trunk
*/

		// Trevor's stuff:

		// mk output
		//VirusTree.printMK();

		// immunity output
		if (Parameters.phenotypeSpace == "geometric") {
			VirusTree.updateRange();
			VirusTree.printRange();
			if (Parameters.immunityReconstruction) {
				printImmunity();
			}
		}		

		// detailed output
		if (Parameters.detailedOutput) {
			printHostPopulation();
		}	

	}



	public void reset() {
		Parameters.day = 0;
		diversity = 0;
		for (int i = 0; i < Parameters.demeCount; i++) {
			HostPopulation hp = demes.get(i);
			hp.reset();
		}
		VirusTree.clear();
	}

}
