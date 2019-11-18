package codaMC;

import java.io.*;
import java.util.*;



public class Merge3RoundsRA {
	
	public static HashMap<String, Double> getMutMap(String file) throws IOException
	{
		ArrayList<String> lines = ReadFile.toList(file);
		HashMap<String, Double> map = new HashMap<String, Double>();
		for(String s : lines) {
			String[] spt = s.trim().split("\t");
			int mutNum = Integer.parseInt(spt[0]);
			if(mutNum == 0) continue;
			if(mutNum > 7) continue;
			ArrayList<Integer> posList = new ArrayList<Integer>();
			for(int i=1;i<mutNum+1;i++) {
				posList.add(Integer.parseInt(spt[i]));
			}
			String mutBaseSeq = spt[mutNum+1];
			if(mutNum == 2) 
				mutBaseSeq += spt[mutNum+2];
			
			String key = "";
			for(int i=0;i<mutNum;i++) {
				key += posList.get(i) + "-";
			}
			key += mutBaseSeq;
			double ra = Double.parseDouble(spt[spt.length-1]);
			map.put(key, ra);
		}
		return map;
	}
	

	public static HashMap<String, Double> getMap(String file) throws IOException {
		String refSeq = "TAACAGGGGGCCACAGCAGAAGCGTTCACGTCGCAGCCCCTGTCAGATTCTGGTGAATCTGCGAATTCTGCTGTATATCTC";
		
		ArrayList<String> lines = ReadFile.toList(file);
		HashMap<String, Double> map = new HashMap<String, Double>();
		for(String s : lines) {
			String[] spt = s.trim().split("\t");
			int mutNum = Integer.parseInt(spt[0]);
			if(mutNum == 0) continue;
			if(mutNum > 7) continue;
			ArrayList<Integer> posList = new ArrayList<Integer>();
			for(int i=1;i<mutNum+1;i++) {
				posList.add(Integer.parseInt(spt[i]));
			}
			String mutBaseSeq = spt[mutNum+1];
			if(mutNum == 2) 
				mutBaseSeq += spt[mutNum+2];
			char[] ref = refSeq.toCharArray();
			for(int k=0;k<mutNum;k++) {
				ref[posList.get(k)] = mutBaseSeq.charAt(k); 
			}
			String seq = new String(ref);
			double ra = Double.parseDouble(spt[spt.length-1]);
			map.put(seq, ra);
		}
		return map;
	}
	
	public static void mergeAll(String file1, String file2, String file3) throws IOException
	{
		HashMap<String, Double> map1 = getMutMap("/export/home/s2982206/rna/scripts/run/C3/var.ra");
		HashMap<String, Double> map2 = getMutMap("/export/home/s2982206/rna/scripts/run/C4/var.ra");
		HashMap<String, Double> map3 = getMutMap("/export/home/s2982206/rna/scripts/run/C5/var.ra");
		
		HashSet<String> keys = new HashSet<String>();
		HashMap<String, Double> map4 = new HashMap<String, Double>();

		for(String s : map1.keySet())
			keys.add(s);
		for(String s : map2.keySet())
			keys.add(s);
		for(String s : map3.keySet())
			keys.add(s);
		for(String key : keys) {
			ArrayList<Double> list = new ArrayList<Double>();
			if(map1.containsKey(key))
				list.add(map1.get(key));
			if(map2.containsKey(key))
				list.add(map2.get(key));
			if(map3.containsKey(key))
				list.add(map3.get(key));
			double mean = math.Stat.mean(list);
			map4.put(key, mean);
		}
		//System.out.println(map4.size());
		for(String key : map4.keySet()) {
			String[] spt = key.split("-");
			if(spt.length == 2)
				System.out.printf("1\t%s\t%s\t-\t100.0\t100.0\t%4.2f\n", spt[0],spt[1],map4.get(key));
		}
		
		for(String key : map4.keySet()) {
			String[] spt = key.split("-");
			if(spt.length == 3) {
				double ra1 = -1;
				double ra2 = -1;
				String key1 = spt[0]+"-"+spt[2].charAt(0);
				String key2 = spt[1]+"-"+spt[2].charAt(1);
				if(map4.containsKey(key1))
					ra1 = map4.get(key1);
				if(map4.containsKey(key2))
					ra2 = map4.get(key2);
				String sra1 = String.format("%4.2f", ra1);
				String sra2 = String.format("%4.2f", ra2);
				if(ra1 < 0) 
					sra1 = "****";
				if(ra2 < 0)
					sra2 = "****";
				System.out.printf("2\t%s\t%s\t%c\t%c\t100.0\t100.0\tnone\t%s\t%s\t%4.2f\n",spt[0],spt[1],spt[2].charAt(0), spt[2].charAt(1),sra1, sra2,map4.get(key));
			}
		}
		
		for(String key : map4.keySet()) {
			String[] spt = key.split("-");
			if(spt.length == 4) {
				double ra1 = -1;
				double ra2 = -1;
				double ra3 = -1;
				String key1 = spt[0]+"-"+spt[3].charAt(0);
				String key2 = spt[1]+"-"+spt[3].charAt(1);
				String key3 = spt[2]+"-"+spt[3].charAt(2);
				if(map4.containsKey(key1))
					ra1 = map4.get(key1);
				if(map4.containsKey(key2))
					ra2 = map4.get(key2);
				if(map4.containsKey(key3))
					ra3 = map4.get(key3);
				String sra1 = String.format("%4.2f", ra1);
				String sra2 = String.format("%4.2f", ra2);
				String sra3 = String.format("%4.2f", ra3);
				if(ra1 < 0) 
					sra1 = "****";
				if(ra2 < 0)
					sra2 = "****";
				if(ra3 < 0)
					sra3 = "****";
				System.out.printf("3\t%s\t%s\t%s\t%s\t100.0\t100.0\t%s\t%s\t%s\tnnn\t****,****,****\t%4.2f\n",spt[0],spt[1],spt[2],spt[3],sra1,sra2,sra3,map4.get(key));
			}
		}
		
		for(String key : map4.keySet()) {
			String[] spt = key.split("-");
			if(spt.length >= 5) {
				double ra = map4.get(key);
				int mutNum = spt.length-1;
				System.out.printf("%d\t", spt.length-1);
				for(int i=0;i<mutNum+1;i++) {
					System.out.printf("%s\t", spt[i]);
				}
				System.out.print("100.0\t100.0\t");
				for(int i=0;i<mutNum;i++) {
					System.out.print("1.00\t");
				}
				for(int i=0;i<mutNum;i++) {
					System.out.print("n");
				}
				System.out.printf("\t********\t%4.2f\n",ra);
			}
		}
		
	}
	
	public static void main(String[] args) throws IOException
	{
		mergeAll(args[0], args[1], args[2]);
	}
}