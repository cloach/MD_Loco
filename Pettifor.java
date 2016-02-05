import java.lang.Math; 
public class Pettifor {
    public static void main(String args[]) {
	    if (args.length != 3) {
		System.out.println("Needs 3 arguments: Name fcc_latice_constnat K_F");
		System.exit(1);
	    }
	    
	    String name = args[0];
	    int N = 4;
	    double fcc_lattice_constant = Double.parseDouble(args[1]);
	    double volume = fcc_lattice_constant*fcc_lattice_constant*fcc_lattice_constant;
	    double K_fermi = Math.pow((3*Math.PI*Math.PI*N/volume),(1.0/3.0));
	    K_fermi = Double.parseDouble(args[2]);
	    //Input Parameters
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
	    //
	    double cutoff = 15.000;
	    int num_samples = 10000;
	    int count = 1;
	    double interval = cutoff/num_samples;
	    
	    System.out.println("#Pot created in java");
	    System.out.println("");
	    System.out.println(name);
	    System.out.println("N " + num_samples);
	    System.out.println("");
	    
	    for (double R = 0.001; R <= cutoff + 0.001; R = R + interval) {
	       	double potential = (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)); 
		double force = (2*Z*Z/(R*R)) * (A_1*Math.cos(K_1*R + a_1)*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (K_1*A_1*Math.sin(K_1*R + a_1)*Math.exp(-kappa_1*R) + K_2*A_2*Math.sin(K_2*R + a_2)*Math.exp(-kappa_2*R) + K_3*A_3*Math.sin(K_3*R + a_3)*Math.exp(-kappa_3*R)) + (2*Z*Z/R) * (A_1*Math.cos(K_1*R + a_1)*kappa_1*Math.exp(-kappa_1*R) + A_2*Math.cos(K_2*R + a_2)*kappa_2*Math.exp(-kappa_2*R) + A_3*Math.cos(K_3*R + a_3)*kappa_3*Math.exp(-kappa_3*R));   
		System.out.println(count + " " + R + " " + potential + " " + force);
		count=count+1;
	    }
    }
}
