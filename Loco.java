import java.io.*; 
import java.util.*;

public class Loco {
    
    public static void main(String args[]) throws IOException {
	//if (args.length != 1) {
	//    System.out.println("Needs 1 filename argument");
	//    System.exit(1);
	//}
	
	double cutoff = 4.0;
	
	double f_e = 1.0;
	double R_0 = 1.0;
	double beta = 1.0;
	double a = 1.0;
	double b = 1.0;
	
	
	String atomList = args[0];
	atomList = atomList.trim();
	//Console myConsole = System.console();
	
	int numLines = countLines(atomList);
	
	//Point latticeA = new Point(2.857,0.0,0.0);
	//Point latticeB = new Point(1.4285,2.474,0.0);
	//Point latticeC = new Point(0.0,0.0,6.997);

	double scale = Double.parseDouble(args[2]);	

	Point latticeA = new Point(1.0,0.0,0.0);
	Point latticeB = new Point(0.0,1.0,0.0);
	Point latticeC = new Point(0.0,0.0,1.0);
	
	//latticeA = latticeA.multiply(scale);
	//latticeB = latticeB.multiply(scale);
	//latticeC = latticeC.multiply(scale);
	
	int supercellA = Integer.parseInt(args[1]);
	int supercellB = Integer.parseInt(args[1]);
	int supercellC = Integer.parseInt(args[1]);
	Point boxA = latticeA.multiply(supercellA);
	Point boxB = latticeB.multiply(supercellB);
	Point boxC = latticeC.multiply(supercellC);
	int numAtoms = numLines * supercellA * supercellB * supercellC;
	
	if ( (Point.magnitude(boxA) <= 2*cutoff) || (Point.magnitude(boxB) < 2*cutoff) || (Point.magnitude(boxA) < 2*cutoff) ) {
	    //System.out.println("A simulation box dimension is not more than 2 times the cutoff! Minimum Images will fuck up!");
	    //System.exit(1);
	}
	
	//fill matrix 
	//double[][] atomMatrix = fillMatrix(atomList, numLines);  
	double[][] atomMatrix = Cell.fillMatrix(atomList, numLines, latticeA, latticeB, latticeC, supercellA, supercellB, supercellC);  
	
	double pairEnergy = 0.0;
	double pairForce = 0.0;
	Point totForce = new Point(0.0,0.0,0.0);
	Point totForceEmbed = new Point(0.0,0.0,0.0);
	//
	double embeddingEnergy = 0.0;
	double embeddingForce = 0.0;
	
	for (int i=0; i<numAtoms; i++) {
	    //System.out.println("i: " + i);
	    Point iForce = new Point(0.0,0.0,0.0);
	    Point iForceEmbed = new Point(0.0,0.0,0.0);
	    double atomPairEnergy = 0.0;
	    double atomPairForce = 0.0;
	    Point iPoint = new Point(atomMatrix[i][0],atomMatrix[i][1],atomMatrix[i][2]);
	    //
	    double chargeDensity = 0.0;
	    for (int j=0; j<numAtoms; j++) {
		if (i != j) {
		    //System.out.println("j: " + j);
		    Point jPoint = new Point(atomMatrix[j][0],atomMatrix[j][1],atomMatrix[j][2]);
		    //double distance = scale*Point.minImage(iPoint, jPoint, boxA, boxB, boxC);
		    
		    //Point vector = Point.minImageVector(iPoint, jPoint, boxA, boxB, boxC);
		    Point vector = Point.subtract(iPoint, jPoint);
		    
		    //System.out.printf("i: %d ",i);
		    //iPoint.write();
		    //System.out.printf("j: %d ",j);
		    //jPoint.write();
		    //System.out.println("R");
		    //vector.write();
		    //System.out.println(" ");
		    double distance = scale*Point.magnitude(vector);
		    
		    //if (distance < cutoff) { 
		    if (1==1) { 
			//double pairPot = pettiforPot(distance);
			double pairPot = lennardJonesPot(distance);
			//double pairPotForce = pettiforForce(distance);
			double pairPotForce = lennardJonesForce(distance);
			Point unitVector = vector.multiply(1/Point.magnitude(vector));
			//System.out.println("should be 1: " + Point.magnitude(unitVector));
			Point jForce = unitVector.multiply(pairPotForce);
			atomPairEnergy = atomPairEnergy + pairPot;
			atomPairForce = atomPairForce + pairPotForce;
			iForce = Point.add(iForce,jForce);
			//System.out.println("jForce: x: " + jForce.getX() + " y: " + jForce.getY() + " z: " + jForce.getZ() + " mag : " + Point.magnitude(jForce));
			//System.out.println("jatomForce: " + pairPotForce);
			//System.out.println("iForce: x: " + iForce.getX() + " y: " + iForce.getY() + " z: " + iForce.getZ() + " mag : " + Point.magnitude(iForce));
			//System.out.println("iatomForce: " + atomPairForce);
			
			//
			double pairCharge = chargeDensity(distance, f_e, R_0, beta);
			double gradRhoMag = gradRho(distance, f_e, R_0, beta);
			Point jForceEmbed = unitVector.multiply(gradRhoMag);
			chargeDensity = chargeDensity + pairCharge;
			iForceEmbed = Point.add(iForceEmbed,jForceEmbed);
		    }
		}
	    }
	    //i

	    
	    
	    pairEnergy = pairEnergy + atomPairEnergy/2.0; 
	    pairForce = pairForce + atomPairForce;
	    totForce = Point.add(totForce,iForce);
	    
	    //update forces for each atom
	    //atomMatrix[i][3] = iForce.getX();
	    //atomMatrix[i][4] = iForce.getY();
	    //atomMatrix[i][5] = iForce.getZ();
	    
	    double atomEmbeddingEnergy = a*chargeDensity - b*Math.pow(chargeDensity,0.5);
	    //double atomEmbeddingForce = (a*beta/R_0)*chargeDensity - (b*beta/(R_0*2.0))*Math.pow(chargeDensity,0.5);
	    double atomEmbeddingForce = (a*beta/R_0)*chargeDensity - (b*beta/(R_0*2.0))*Math.pow(chargeDensity,0.5);
	    atomMatrix[i][3] = atomEmbeddingForce;
	    atomMatrix[i][4] = -chargeDensity/R_0;

	    embeddingEnergy = embeddingEnergy + atomEmbeddingEnergy;
	    embeddingForce = embeddingForce + atomEmbeddingForce;
	    
	    double dFdRho = (a-1/2*beta/Math.pow(chargeDensity,0.5)) ;
	    Point atomEmbeddingForceVector = iForceEmbed.multiply(dFdRho);
	    
	    /////REOVE THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    //atomMatrix[i][3] = atomEmbeddingForceVector.getX();
	    //atomMatrix[i][4] = atomEmbeddingForceVector.getY();
	    //atomMatrix[i][5] = atomEmbeddingForceVector.getZ();
	    ////////////////////////////////////////////
	    
	    totForceEmbed = Point.add(totForceEmbed,atomEmbeddingForceVector); 
	    
	    System.out.println("nonvectorforce: " + atomEmbeddingForce);
	    System.out.println("nonvectorforce: " + Point.magnitude(atomEmbeddingForceVector));
	}
	
	Point totForceGlue = new Point(0.0,0.0,0.0);	
	for (int i=0; i<numAtoms; i++) {
	    Point iPoint = new Point(atomMatrix[i][0],atomMatrix[i][1],atomMatrix[i][2]);
	    Point iForceGlue = new Point(0.0,0.0,0.0);
	    //double chargeDensity = 0.0;
	    for (int j=0; j<numAtoms; j++) {
		if (i!=j) {
		    Point jPoint = new Point(atomMatrix[j][0],atomMatrix[j][1],atomMatrix[j][2]);
		    Point vector = Point.subtract(iPoint, jPoint);
		    double distance = scale*Point.magnitude(vector);
		    Point unitVector = vector.multiply(1/Point.magnitude(vector));
		    double glueForcejMag = -( (atomMatrix[i][3]*atomMatrix[j][4] + atomMatrix[j][3]*atomMatrix[i][4]) );
		    Point jForceGlue = unitVector.multiply(glueForcejMag);
		    //double distance = scale*Point.minImage(iPoint, jPoint, boxA, boxB, boxC);
	            //double pairCharge = chargeDensity(distance, f_e, R_0, beta);
		    iForceGlue = Point.add(iForceGlue,jForceGlue);
		}
	    }
	    totForceGlue = Point.add(totForceGlue,iForceGlue);
	    System.out.println("magx" + i + ": " + +iForceGlue.getX());
	}
	
	
	
	
       	Cell.writeMatrix(atomMatrix);
	System.out.println(" ");
	
	double pairEnergyPerAtom = pairEnergy/numAtoms;    
	double pairForcePerAtom = pairForce/numAtoms;    
	System.out.println("pairEnergy = " + pairEnergy);
	System.out.println("pairEnergyPerAtom = " + pairEnergyPerAtom);
	System.out.println("pairForce = " + pairForce);
	System.out.println("pairForcePerAtom = " + pairForcePerAtom);
	System.out.println("totPairForceVector: x: " + totForce.getX() + " y: " + totForce.getY() + " z: " + totForce.getZ() + " mag: " + Point.magnitude(totForce));
	System.out.println("totPairForceVectorMag = " + Point.magnitude(totForce));
	System.out.println("totalEmbedVectorMag = " + Point.magnitude(totForceEmbed));
	
	System.out.printf("totalEmbedVector: ");
	totForceEmbed.write();
	
	//
	double embeddingEnergyPerAtom = embeddingEnergy/numAtoms;
	double embeddingForcePerAtom = embeddingForce/numAtoms;
	System.out.println("EmbeddingEnergy = " + embeddingEnergy);
	System.out.println("EmbeddingEnergyPerAtom = " + embeddingEnergyPerAtom);
	System.out.println("EmbeddingForce = " + embeddingForce);
	System.out.println("EmbeddingForcePerAtom = " + embeddingForcePerAtom);
	double totalEnergy = pairEnergy + embeddingEnergy;
	double totalForce = pairForce + embeddingForce;
	double totalEnergyPerAtom = totalEnergy/numAtoms;
	double totalForcePerAtom = totalForce/numAtoms;
	System.out.println("TotalEnergy = " + totalEnergy);
	System.out.println("TotalEnergyPerAtom = " + totalEnergyPerAtom);
	System.out.println("TotalForce = " + totalForce);
	System.out.println("TotalForcePerAtom = " + totalForcePerAtom);
	
	System.exit(0);	
    }
    
    public static int countLines(String filename) throws IOException {
	InputStream is = new BufferedInputStream(new FileInputStream(filename));
	try {
	    byte[] c = new byte[1024];
	    int count = 0;
	    int readChars = 0;
	    boolean empty = true;
	    while ((readChars = is.read(c)) != -1) {
		empty = false;
		for (int i = 0; i < readChars; ++i) {
		    if (c[i] == '\n') {
			++count;
		    }
		}
	    }
	    return (count == 0 && !empty) ? 1 : count;
	} finally {
	    is.close();
	}
    }
    
    public static double[][] fillMatrix(String filename, int matrixSize) throws IOException {
	double[][] atomMatrix = new double[matrixSize][3]; 
	BufferedReader atomBuffer = new BufferedReader(new FileReader(filename));
	Scanner scan = new Scanner(atomBuffer);
	int countLine = 1;
	while(scan.hasNext()) {
	    String atomLine = scan.nextLine();
	    atomLine = atomLine.trim();
	    String[] atomLineList = atomLine.split(" ");
	    int atomNumber = Integer.parseInt(atomLineList[0]);
	    double xValue = Double.parseDouble(atomLineList[1]);
	    double yValue = Double.parseDouble(atomLineList[2]);
	    double zValue = Double.parseDouble(atomLineList[3]);
	    atomMatrix[countLine][0] = xValue;
	    atomMatrix[countLine][1] = yValue;
	    atomMatrix[countLine][2] = zValue;
	    countLine++;
	}
	return atomMatrix;
    }
    
    public static double pettiforPot(double R) {
	double K_fermi = 1.8;
	double Z = 1.0;
	double A_1 = 7.964;
	double A_2 = 1.275;
	double A_3 = 0.030;
	double K_1 = 2*K_fermi * 0.156;
	double K_2 = 2*K_fermi * 0.644;
	double K_3 = 2*K_fermi * 0.958;
	double kappa_1 = 2*K_fermi * 0.793;
	double kappa_2 = 2*K_fermi * 0.698;
	double kappa_3 = 2*K_fermi * 0.279;
	double a_1 = Math.PI * -0.441;
	double a_2 = Math.PI * 0.832;
	double a_3 = Math.PI * 0.431;
	
	double potential = (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)); 
	double force = (2*Z*Z/(R*R)) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (K_1*A_1*Math.sin(K_1*R + a_1)*Math.exp(-kappa_1*R) + K_2*A_2*Math.sin(K_2*R + a_2)*Math.exp(-kappa_2*R) + K_3*A_3*Math.sin(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*kappa_1*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*kappa_2*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*kappa_3*Math.exp(-kappa_3*R));   

	return potential;
    }

    public static double pettiforForce(double R) {
	double K_fermi = 1.8;
	double Z = 1.0;
	double A_1 = 7.964;
	double A_2 = 1.275;
	double A_3 = 0.030;
	double K_1 = 2*K_fermi * 0.156;
	double K_2 = 2*K_fermi * 0.644;
	double K_3 = 2*K_fermi * 0.958;
	double kappa_1 = 2*K_fermi * 0.793;
	double kappa_2 = 2*K_fermi * 0.698;
	double kappa_3 = 2*K_fermi * 0.279;
	double a_1 = Math.PI * -0.441;
	double a_2 = Math.PI * 0.832;
	double a_3 = Math.PI * 0.431;
	
	double potential = (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)); 
	double force = (2*Z*Z/(R*R)) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (K_1*A_1*Math.sin(K_1*R + a_1)*Math.exp(-kappa_1*R) + K_2*A_2*Math.sin(K_2*R + a_2)*Math.exp(-kappa_2*R) + K_3*A_3*Math.sin(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*kappa_1*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*kappa_2*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*kappa_3*Math.exp(-kappa_3*R));   
	
	return force;
    }
    
    public static double chargeDensity(double R, double f_e, double R_0, double beta) {
	double potential = f_e*Math.exp(-beta*(R/R_0 - 1)); 	
	
	return potential;
    }

    public static double gradRho(double R, double f_e, double R_0, double beta) {
	double potential = -f_e*beta/R_0*Math.exp(-beta*(R/R_0 - 1)); 	
	
	return potential;
    }
    
    public static double lennardJonesPot(double R) {
	double R_0 = 0.707; //change this later
	double potential = Math.pow(R_0/R,12) - 2*Math.pow(R_0/R,6);
	return potential;
    }

    public static double lennardJonesForce(double R) {
	double R_0 = 0.707; //change tihs later
	double force = (12/R_0)*(Math.pow(R_0/R,13) - Math.pow(R_0/R,7));
	return force;
    }
				     
}


class simulation {
    double[][] atomMatrix;
    double density, potenetialEnergy, kineticEnergy
