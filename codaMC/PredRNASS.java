package codaMC;
import java.io.*;
import java.util.*;

public class PredRNASS {

	private int stepNum = 500000;
	private double T0 = 10.0;
	private double T1 = 0.1;
	private String seq;
	private RNAEnergyCalculator ec;
	
	public PredRNASS(String seq) {
		this.seq = seq;
		ec = new RNAEnergyCalculator(seq);
	}

        public static ArrayList<String> toList(String fileName) throws IOException
        {
                BufferedReader in = new BufferedReader(new FileReader(fileName));
                ArrayList<String> list = new ArrayList<String>();
                String line;
                while((line = in.readLine())!=null)
                {
                        if(in.ready() || line.length() > 0)
                                list.add(line);
                }
                in.close();

                return list;
        }

	
	public PredRNASS(String seq, String mtxFile, double lamda) throws IOException {
		this.seq = seq;
		ec = new RNAEnergyCalculator(seq);
		ec.addXpPredInfo(mtxFile, lamda);
	}
	
	public boolean accept(double mutEnergy, double T, Random rand)
	{
		if(T == 0)
		{
			if(mutEnergy >0)
				return false;
			else
				return true;
		}
		else if(mutEnergy > 0)
		{
			double pAc = Math.exp(-1.0*mutEnergy/T);
			return rand.nextDouble() < pAc;
		}
		else
			return true;
	}
	
	
	
	public RNASSUnit mc() {
		

		RNASSUnit unit = new RNASSUnit(seq.length());
		double lastE = ec.getEnergy(unit);
		
		RNASSUnit bestResult = unit.clone();
		double bestE = ec.getEnergy(unit);
		Random rand = new Random();
		
		for(double T = T0;T > T1;T=T*0.9) {
			for(int k=0;k<stepNum;k++) {
				RNASSUnit curUnit = unit.randomMove(rand);
				double curE = ec.getEnergy(curUnit);
				if(accept(curE-lastE, T, rand)) {
					lastE = curE;
					unit = curUnit;
					if(curE < bestE) {
						bestResult = curUnit.clone();
						bestE = curE;
					}
				}
			}
		}
		//System.out.println(bestResult.toSecSequence() + " " + bestE);
		
		
		return bestResult;
	}
	

	
	public static void main(String[] args) throws IOException
	{
		String seq = args[0];
		
		
		PredRNASS ss = new PredRNASS(seq.toUpperCase());

		RNASSUnit result = ss.mc();
		String sec = result.toSecSequence();
		System.out.println(sec);
	}

}
