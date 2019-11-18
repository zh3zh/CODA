package codaMC;
import java.io.*;
import java.util.*;

public class RNASSUnit {

	public int len;
	public int[] pairIndexList;

	
	public RNASSUnit(int len) {
		this.len = len;
		this.pairIndexList = new int[len];
		for(int i=0;i<len;i++) {
			this.pairIndexList[i] = -1;
		}
	}
	
	public RNASSUnit(int len, int[] indexList) {
		this.len = len;
		this.pairIndexList = new int[len];
		for(int i=0;i<len;i++) {
			this.pairIndexList[i] = indexList[i];
		}
	}
	
	public void addPair(int posA, int posB) {
		if(this.pairIndexList[posA] >= 0) {
			int pA = pairIndexList[posA];
			this.pairIndexList[posA] = -1;
			this.pairIndexList[pA] = -1;
		}
		else if(this.pairIndexList[posB] >= 0) {
			int pB = pairIndexList[posB];
			this.pairIndexList[posB] = -1;
			this.pairIndexList[pB] = -1;
		}
		this.pairIndexList[posA] = posB;
		this.pairIndexList[posB] = posA;
	}
	
	public RNASSUnit clone() {
		return new RNASSUnit(len, pairIndexList);
	}
	
	public void printMtx(String file) throws IOException {
		PrintStream out = new PrintStream(file);
		for(int i=0;i<len;i++) {
			for(int j=0;j<len-1;j++) {
				double s = 0.0;
				if(pairIndexList[i] == j)
					s = 1.0;
				out.printf("%5.3f\t", s);
			}
			double s = 0.0;
			if(pairIndexList[i] == len-1)
				s = 1.0;
			out.printf("%5.3f\n", s);
		}
		out.close();
	}
	
	public void printContact(String seq) {
		for(int i=0;i<len;i++) {
			if(pairIndexList[i] < 0) continue;
			int pi = pairIndexList[i];
			System.out.printf("%c%d-%c%d\n",seq.charAt(i), i, seq.charAt(pi), pi);
		}
	}
	
    public String toSecSequence() {
    	HashMap<Character, Integer> map = new HashMap<Character, Integer>();
    	map.put('(', -1);
    	map.put(')', 1);
    	map.put('[', -100);
    	map.put(']', 100);
    	map.put('{', -10000);
    	map.put('}', 10000);
    	map.put('<', -1000000);
    	map.put('>', 1000000);
    	map.put('.', 0);
    	ArrayList<String> bractList = new ArrayList<String>();
    	bractList.add("()");
    	bractList.add("{}");
    	bractList.add("[]");
    	
    	char[] ss = new char[len];
    	for(int i=0;i<len;i++) {
    		ss[i] = '.';
    	}
    	for(int i=0;i<len;i++) {
    		int pi = pairIndexList[i];
    		if(pi < i) continue;
    		if(pi == -1) continue;
    		int idA = i;
    		int idB = pi;
    		int x = 0;
    		int y = 0;
    		int z = 0;
    		for(int k=idA+1;k<idB;k++) {
    			char c = ss[k];
    			if(c == '(') x++;
    			if(c == ')') x--;
    			if(c == '[') y++;
    			if(c == ']') y--;
    			if(c == '{') z++;
    			if(c == '}') z--;
    		}
    		if(x == 0) {
    			ss[idA] = '(';
    			ss[idB] = ')';
    		}
    		else if(y == 0) {
    			ss[idA] = '[';
    			ss[idB] = ']';
    		}
    		else if(z == 0){
    			ss[idA] = '{';
    			ss[idB] = '}';
    		}
    		else {
    			ss[idA] = '<';
    			ss[idB] = '>';
    		}
    	}
    	return new String(ss);
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

	
	public RNASSUnit randomMove(Random rand) {
		
		int[] newPairIndexList = new int[len];
		for(int i=0;i<len;i++) {
			newPairIndexList[i] = pairIndexList[i];
		}
		
		int posA = rand.nextInt(len);
		int pA = pairIndexList[posA];
		if(pA == -1) {
			int posB = rand.nextInt(len);
			while(posB == posA)
				posB = rand.nextInt(len);
			int pB = pairIndexList[posB];
			
			if(pB == -1) {
				newPairIndexList[posA] = posB;
				newPairIndexList[posB] = posA;
			}
			else if(rand.nextBoolean()) {
				newPairIndexList[pB] = -1;
				newPairIndexList[posA] = posB;
				newPairIndexList[posB] = posA;
			}
		}
		else {
			int posB = pA;
			if(rand.nextBoolean()) { 
				newPairIndexList[posA] = -1;
				newPairIndexList[posB] = -1;
			}
		}
		return new RNASSUnit(len, newPairIndexList);
	}
	
	public double accuracy(String trueContactMtx) throws IOException {
		ArrayList<String> lines = toList(trueContactMtx);
		double[][] contact = new double[len][len];
		for(int i=0;i<len;i++) {
			String s = lines.get(i);
			String[] spt = s.trim().split("\t");
			for(int j=0;j<len;j++) {
				contact[i][j] = Double.parseDouble(spt[j]);
			}
		}
		
		int tp = 0;
		int fp = 0;
		int tn = 0;
		int fn = 0;

		for(int i=0;i<len;i++) {
			for(int j=i+1;j<len;j++) {

				if(contact[i][j] > 0 && pairIndexList[i] == j)
					tp++;
				else if(contact[i][j] > 0 && pairIndexList[i] != j)
					fn++;
				else if(pairIndexList[i] == j)
					fp++;
				else if(pairIndexList[i] != j)
					tn++;
			}
		}
		
		double N = tn + tp + fn + fp;
		double S = (tp+fn)*1.0/N;
		double P = (tp+fp)*1.0/N;
		double MCC = (tp/N - S*P)/Math.sqrt(P*S*(1-P)*(1-S));
		return MCC;
	}
	
}
