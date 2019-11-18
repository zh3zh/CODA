package codaMC;
import java.io.*;
import java.util.*;


public class PredictSSToMtx {
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


	public static void readSS(String input, String output) throws IOException
	{
		int len = 81;
		double[][] mtx = new double[81][81];
		double fixValue = 10.0;
		mtx[5][40] = fixValue;
		mtx[6][39] = fixValue;
		mtx[7][38] = fixValue;
		mtx[8][37] = fixValue;
		mtx[9][36] = fixValue;
		mtx[10][35] = fixValue;
		mtx[11][34] = fixValue;
		mtx[14][71] = fixValue;
		mtx[15][70] = fixValue;
		mtx[16][69] = fixValue;
		mtx[17][68] = fixValue;
		mtx[18][67] = fixValue;
		mtx[19][66] = fixValue;
		mtx[20][65] = fixValue;
		mtx[21][33] = fixValue;
		mtx[22][32] = fixValue;
		mtx[23][31] = fixValue;
		mtx[25][42] = fixValue;
		mtx[26][41] = fixValue;
		mtx[43][60] = fixValue;
		mtx[44][59] = fixValue;
		mtx[45][58] = fixValue;
		mtx[46][57] = fixValue;
		mtx[47][56] = fixValue;
		mtx[48][55] = fixValue;
		mtx[49][54] = fixValue;
		
		ArrayList<String> ssList= toList(input);
		for(String s : ssList) {
			for(int i=0;i<s.length();i++) {
				if(s.charAt(i) == '(') {
					int x = 0;
					for(int j=i+1;j<s.length();j++) {
						char c = s.charAt(j);
						if(c == '(')
							x ++;
						if(c == ')') {
							if(x == 0) {
								mtx[j][i] += 0.1;
								break;
							}
							else
								x--;
						}
					}
				}
				
				if(s.charAt(i) == '[') {
					int x = 0;
					for(int j=i+1;j<s.length();j++) {
						char c = s.charAt(j);
						if(c == '[')
							x ++;
						if(c == ']') {
							if(x == 0) {
								mtx[j][i] += 0.1;
								break;
							}
							else
								x--;
						}
					}
				}
				
				if(s.charAt(i) == '{') {
					int x = 0;
					for(int j=i+1;j<s.length();j++) {
						char c = s.charAt(j);
						if(c == '{') 
							x ++;
						if(c == '}') {
							if(x == 0) {
								mtx[j][i] += 0.1;
								break;
							}
							else
								x--;
						}
					}
				}
				
				if(s.charAt(i) == '<') {
					int x = 0;
					for(int j=i+1;j<s.length();j++) {
						char c = s.charAt(j);
						if(c == '<')
							x ++;
						if(c == '>') {
							if(x == 0) {
								mtx[j][i] += 0.1;
								break;
							}
							else
								x--;
						}
					}
				}
				
				
			}
		}
		PrintStream out = new PrintStream(output);
		for(int i=0;i<len;i++) {
			for(int j=0;j<len-1;j++) {
				out.printf("%-8.6f\t", mtx[i][j]);
			}
			out.printf("%8.6f\n",mtx[i][len-1]);
		}
		out.close();
		
	}
	
	public static void main(String[] args) throws IOException{
		readSS(args[0], args[1]);
	}
}
