package codaMC;

import java.io.*;
import java.util.*;

public class RunPredict {
	public static void main(String[] args) throws IOException
	{
		String seq = args[0];
		String mtxFile = args[1];
		double wtFactor = Double.parseDouble(args[2]);
		PredRNASS ps = new PredRNASS(seq, mtxFile, wtFactor);
		
		HashMap<String, Integer> ssMap = new HashMap<String, Integer>();
		for(int i=0;i<100;i++) {
			RNASSUnit su = ps.mc();
			String ss = su.toSecSequence();
			System.out.println(ss);
			if(ssMap.containsKey(ss))
				ssMap.put(ss, ssMap.get(ss)+1);
			else
				ssMap.put(ss, 1);
		}
	}
}
