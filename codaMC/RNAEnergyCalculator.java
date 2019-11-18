package codaMC;
import java.io.*;
import java.util.*;

public class RNAEnergyCalculator {
	private int len;
	private String seq;
	private HashMap<String, Double> energyMap;
	private HashSet<String> pairSet;
	public double[][] scoreMatrix;
	
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


	public RNAEnergyCalculator(String seq) {
		
		char[] s = seq.toCharArray();
		for(int i=0;i<s.length;i++) {
			if(s[i] == 'T')
				s[i] = 'U';
		}
		
		this.seq = new String(s);
		
		this.len = seq.length();
		this.scoreMatrix = new double[len][len];

		this.energyMap = new HashMap<String, Double>();
		this.energyMap.put("AAUU", -0.93);
		this.energyMap.put("UUAA", -0.93);
		this.energyMap.put("AUUA", -1.1);
		this.energyMap.put("UAAU", -1.33);
		this.energyMap.put("CUGA", -2.08);
		this.energyMap.put("AGUC", -2.08);
		this.energyMap.put("CAGU", -2.21);
		this.energyMap.put("UGAC", -2.21);
		this.energyMap.put("GUCA", -2.24);
		this.energyMap.put("ACUG", -2.24);
		this.energyMap.put("GACU", -2.35);
		this.energyMap.put("UCAG", -2.35);
		this.energyMap.put("CGGC", -2.36);
		this.energyMap.put("GGCC", -3.26);
		this.energyMap.put("CCGG", -3.26);
		this.energyMap.put("GCCG", -3.42);
		
		this.pairSet = new HashSet<String>();
		pairSet.add("AU");
		pairSet.add("UA");
		pairSet.add("GC");
		pairSet.add("CG");
		pairSet.add("GU");
		pairSet.add("UG");
		
	}
	
	public RNAEnergyCalculator(String seq, String dcaFile, double dcaWeight) throws IOException {
		this.seq = seq;
		this.len = seq.length();
		this.scoreMatrix = new double[len][len];

		this.energyMap = new HashMap<String, Double>();
		this.energyMap.put("AAUU", -0.93);
		this.energyMap.put("UUAA", -0.93);
		this.energyMap.put("AUUA", -1.1);
		this.energyMap.put("UAAU", -1.33);
		this.energyMap.put("CUGA", -2.08);
		this.energyMap.put("AGUC", -2.08);
		this.energyMap.put("CAGU", -2.21);
		this.energyMap.put("UGAC", -2.21);
		this.energyMap.put("GUCA", -2.24);
		this.energyMap.put("ACUG", -2.24);
		this.energyMap.put("GACU", -2.35);
		this.energyMap.put("UCAG", -2.35);
		this.energyMap.put("CGGC", -2.36);
		this.energyMap.put("GGCC", -3.26);
		this.energyMap.put("CCGG", -3.26);
		this.energyMap.put("GCCG", -3.42);
		
		this.pairSet = new HashSet<String>();
		pairSet.add("AU");
		pairSet.add("UA");
		pairSet.add("GC");
		pairSet.add("CG");
		pairSet.add("GU");
		pairSet.add("UG");
		
		ArrayList<String> lines = toList(dcaFile);
		for(String line : lines) {
			String[] spt = line.trim().split(",");
			int idA = Integer.parseInt(spt[0])-1;
			int idB = Integer.parseInt(spt[1])-1;
			double score = Double.parseDouble(spt[2]);
			scoreMatrix[idA][idB] = -dcaWeight * score;
		}
	}
	
	
	public void addXpPredInfo(String predMtx, double lamda) throws IOException {
		
		ArrayList<Double> scores = new ArrayList<Double>();
		ArrayList<Integer> IList = new ArrayList<Integer>();
		ArrayList<Integer> JList = new ArrayList<Integer>();
		ArrayList<String> lines = toList(predMtx);
		int len = 0;
		for(int i=0;i<lines.size();i++) {
			String line = lines.get(i);
			String[] spt = line.trim().split("\t");
			if(len == 0)
				len = spt.length;
			for(int j=0;j<i-3;j++) {
				scores.add(Double.parseDouble(spt[j]));
				IList.add(i);
				JList.add(j);
			}
		}

		double tot = 0;
		for(int i=0;i<scores.size();i++) {
			tot += scores.get(i)*scores.get(i);
		}
		
		double mean = Math.sqrt(tot/scores.size());
		double wt = lamda/mean;
		System.out.println("matrix weight: " + wt);
		for(int i=0;i<scores.size();i++) {
			scoreMatrix[IList.get(i)][JList.get(i)] = -wt * scores.get(i);
			scoreMatrix[JList.get(i)][IList.get(i)] = -wt * scores.get(i);
		}
	}

	public int helixLengthLeft(RNASSUnit unit, int index) {
		int n1 = 0;
		if(index>0) {
			int piLeft = unit.pairIndexList[index-1];
			if(piLeft >=0) {
				n1 = 1;
				for(int i=index-2;i>=0;i--) {
					if(unit.pairIndexList[i] == piLeft + index-1-i)
						n1++;
					else
						break;
				}
			}
		}
		return n1;
	}
	
	public int helixLengthRight(RNASSUnit unit, int index) {
		
		int n2 = 0;
		if(index < len-1) {
			int piRight = unit.pairIndexList[index+1];
			if(piRight >=0) {
				n2 = 1;
				for(int i=index+2;i<len;i++) {
					if(unit.pairIndexList[i] == piRight + i - index-1)
						n2++;
					else
						break;
				}
			}
		}
		return n2;
	}
	
	public double getEnergy(RNASSUnit unit) {
		double e = 0.0;
		for(int i=0;i<len;i++) {
			int pi = unit.pairIndexList[i];
			if(pi == -1) continue;
			if(pi <= i) continue;
			int sep = Math.abs(pi-i);
			if(sep < 4)
				e += 100.0;
			
			String pair = "" + seq.charAt(i) + seq.charAt(pi);
			
			if(!pairSet.contains(pair)) {
				e += 4.0;
				if(i>0 && unit.pairIndexList[i-1] == pi+1 && helixLengthLeft(unit, i) > 3)
					e -= 3.0;
				if(i<len-1 && unit.pairIndexList[i+1] == pi -1 && helixLengthRight(unit, i) > 3)
					e -= 3.0;
			}
				
			
			e += scoreMatrix[i][pi];
			
			/*
			 * single pair penalty: 4
			 * AU and penalty: 0.45
			 */
			
			double penalty = 4.0;
			double auPenalty = 0.45;
			if(i == 0) {
				if(unit.pairIndexList[i+1] != pi-1)
				{
					e += penalty;
					if(pair.equals("AU") || pair.equals("UA")) 
						e += auPenalty;
				}
			}
			else if(i==len-1) {
				if(unit.pairIndexList[i-1] != pi+1)
				{
					e += penalty;
					if(pair.equals("AU") || pair.equals("UA"))
						e += auPenalty;
				}
			}
			else if(unit.pairIndexList[i-1] != pi+1 && unit.pairIndexList[i+1] != pi-1)
			{
				e += penalty;
				if(pair.equals("AU") || pair.equals("UA"))
					e += auPenalty*2;
			}
			else if(unit.pairIndexList[i-1] != pi+1)
			{
				if(pair.equals("AU") || pair.equals("UA"))
					e += auPenalty;
			}
			else if(unit.pairIndexList[i+1] != pi-1)
			{
				if(pair.equals("AU") || pair.equals("UA"))
					e += auPenalty;
			}
			
			if(i == 0) continue;
			
			int pj = unit.pairIndexList[i-1];
			
			if(pj == -1) {
				continue;
			}
			
			String key = "" + seq.charAt(i-1)+seq.charAt(i)+seq.charAt(pj)+seq.charAt(pi);
			if(pj != pi+1) {
				continue;
			}
			if(energyMap.containsKey(key))
				e += energyMap.get(key);
		}
		return e;
	}
	
	public static void main(String[] args) throws IOException
	{
		RNASSUnit unit = new RNASSUnit(81);
		unit.pairIndexList[0] = -1;
		unit.pairIndexList[1] = -1;
		unit.pairIndexList[2] = -1;
		unit.pairIndexList[3] = -1;
		unit.pairIndexList[4] = -1;
		unit.pairIndexList[5] = -1;
		unit.pairIndexList[6] = -1;
		unit.pairIndexList[7] = 38;
		unit.pairIndexList[8] = 37;
		unit.pairIndexList[9] = 36;
		unit.pairIndexList[10] = 35;
		unit.pairIndexList[11] = -1;
		unit.pairIndexList[12] = -1;
		unit.pairIndexList[13] = -1;
		unit.pairIndexList[14] = 71;
		unit.pairIndexList[15] = 70;
		unit.pairIndexList[16] = 69;
		unit.pairIndexList[17] = 68;
		unit.pairIndexList[18] = 67;
		unit.pairIndexList[19] = 66;
		unit.pairIndexList[20] = 65;
		unit.pairIndexList[21] = -1;
		unit.pairIndexList[22] = 32;
		unit.pairIndexList[23] = 31;
		unit.pairIndexList[24] = -1;
		unit.pairIndexList[25] = 42;
		unit.pairIndexList[26] = 41;
		unit.pairIndexList[27] = -1;
		unit.pairIndexList[28] = -1;
		unit.pairIndexList[29] = -1;
		unit.pairIndexList[30] = -1;
		unit.pairIndexList[31] = 23;
		unit.pairIndexList[32] = 22;
		unit.pairIndexList[33] = -1;
		unit.pairIndexList[34] = -1;
		unit.pairIndexList[35] = 10;
		unit.pairIndexList[36] = 9;
		unit.pairIndexList[37] = 8;
		unit.pairIndexList[38] = 7;
		unit.pairIndexList[39] = -1;
		unit.pairIndexList[40] = -1;
		unit.pairIndexList[41] = 26;
		unit.pairIndexList[42] = 25;
		unit.pairIndexList[43] = 60;
		unit.pairIndexList[44] = 59;
		unit.pairIndexList[45] = 58;
		unit.pairIndexList[46] = 57;
		unit.pairIndexList[47] = 56;
		unit.pairIndexList[48] = 55;
		unit.pairIndexList[49] = 54;
		unit.pairIndexList[50] = -1;
		unit.pairIndexList[51] = -1;
		unit.pairIndexList[52] = -1;
		unit.pairIndexList[53] = -1;
		unit.pairIndexList[54] = 49;
		unit.pairIndexList[55] = 48;
		unit.pairIndexList[56] = 47;
		unit.pairIndexList[57] = 46;
		unit.pairIndexList[58] = 45;
		unit.pairIndexList[59] = 44;
		unit.pairIndexList[60] = 43;
		unit.pairIndexList[61] = -1;
		unit.pairIndexList[62] = -1;
		unit.pairIndexList[63] = -1;
		unit.pairIndexList[64] = -1;
		unit.pairIndexList[65] = 20;
		unit.pairIndexList[66] = 19;
		unit.pairIndexList[67] = 18;
		unit.pairIndexList[68] = 17;
		unit.pairIndexList[69] = 16;
		unit.pairIndexList[70] = 15;
		unit.pairIndexList[71] = 14;
		unit.pairIndexList[72] = -1;
		unit.pairIndexList[73] = -1;
		unit.pairIndexList[74] = -1;
		unit.pairIndexList[75] = -1;
		unit.pairIndexList[76] = -1;
		unit.pairIndexList[77] = -1;
		unit.pairIndexList[78] = -1;
		unit.pairIndexList[79] = -1;
		unit.pairIndexList[80] = -1;
		
		String seq = "UAACAGGGGGCCACAGCAGAAGCGUUCACGUCGCAGCCCCUGUCAGAUUCUGGUGAAUCUGCGAAUUCUGCUGUAUAUCUC";
		String mtxFile = "C:\\work\\cpeb\\pred.mtx";
		RNAEnergyCalculator ec = new RNAEnergyCalculator(seq);
		ec.addXpPredInfo(mtxFile, 1.0);
		double e = ec.getEnergy(unit);
		System.out.println(e);
		
		unit.pairIndexList[5] = 40;
		unit.pairIndexList[6] = 39;
		unit.pairIndexList[7] = 38;
		unit.pairIndexList[8] = 37;
		unit.pairIndexList[9] = 36;
		unit.pairIndexList[10] = 35;
//		unit.pairIndexList[11] = 34;
		
		unit.pairIndexList[40] = 5;
		unit.pairIndexList[39] = 6;
		unit.pairIndexList[38] = 7;
		unit.pairIndexList[37] = 8;
		unit.pairIndexList[36] = 9;
		unit.pairIndexList[35] = 10;
//		unit.pairIndexList[34] = 11;
		
		double x1 = ec.scoreMatrix[5][40];
		double x2 = ec.scoreMatrix[6][39];
		double x3 = ec.scoreMatrix[7][38];
		double x4 = ec.scoreMatrix[8][37];
		double x5 = ec.scoreMatrix[9][36];
		double x6 = ec.scoreMatrix[10][35];
		double x7 = ec.scoreMatrix[11][34];
		

		System.out.println(x1);
		System.out.println(x2);
		System.out.println(x3);
		System.out.println(x4);
		System.out.println(x5);
		System.out.println(x6);
		System.out.println(x7);
		
		
		e = ec.getEnergy(unit);
		System.out.println(e);
		
		
	}
}
