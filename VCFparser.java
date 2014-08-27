//--- Java imports ---
import java.util.*;
import java.util.zip.*;
import java.io.*;
import java.lang.reflect.*;

class VCFparser
{
    private FileInputStream _fin   = null;
    private BufferedReader _br     = null;
    private String         _line   = "";
    private int            _id_ind = -1;
    private String _id = "";
    private boolean _pass = false;
    
    public VCFparser(String fileName,String id,boolean pass)
    {
	try {
	    _fin = new FileInputStream(fileName);
	    InputStream is = _fin;
	    String tmp = new String(fileName);
	    tmp.toLowerCase();
	    if (tmp.endsWith(".gz") || tmp.endsWith(".gzip"))
		is = new GZIPInputStream(_fin);
 	    _br = new BufferedReader(new InputStreamReader(is));
	    readLine();
	} catch (Exception ex) {
	    System.err.println("Can't open file " + fileName);
	    System.err.println(ex.toString());
	}
	if (id != null) _id = id;
	_pass = pass;
    }
    
    public boolean hasMoreInput()
    {
	return (_line != null);
    }
    
    public Variant parseLine()
    {
	String line = _line;
	readLine();
	
	if (line == null || line.length() == 0) return null;

	StringTokenizer toks = new StringTokenizer(line);
	if (line.startsWith("#")) {
	    if (line.startsWith("#CHROM")) {
		int index = 0;
		while (toks.hasMoreTokens()) {
		    index++;
		    String tok = toks.nextToken();
		    if (tok.equals(_id)) _id_ind = index;
		}
	    }
	    return null;
	}


	int index = 0,genotype_ind = -1;
	int chr = -1,pos = -1;
	String REF = "",FILTER = "",ALT = "",phase = "0/0",INFO = "";
	while (toks.hasMoreTokens()) {
	    index++;
	    if (index == 1) // Parsing chromosome
		chr = getChromIndex(toks.nextToken());
	    else if (index == 2) // Parsing position
		pos = Integer.parseInt(toks.nextToken());
	    else if (index == 4) // Parsing reference allele
		REF = toks.nextToken();
	    else if (index == 5) // Parsing alternative allele
		ALT = toks.nextToken();
	    else if (index == 7) // FILTER field
		FILTER = toks.nextToken();
// 	    else if (index == 8) // INFO field
// 		INFO = toks.nextToken();
	    else if (index == 9) { // Output format
		genotype_ind = getGenotypeIndex(toks.nextToken());
		if (genotype_ind < 0) break;
	    } else if (index == _id_ind) { // Phasing
		StringTokenizer words =
		    new StringTokenizer(toks.nextToken(),":");
		phase = words.nextToken();
		break;
	    } else toks.nextToken();
	}

	if (ALT.equals("<DEL>"))                 return null; // Imprecise SV
	if (_pass && FILTER.indexOf("PASS") < 0) return null; // Filtered out

	// Upper casing
	REF = REF.toUpperCase();
	ALT = ALT.toUpperCase();    

	// Splitting
	String[] alts = ALT.split(",");
	int n = alts.length;
	
	// Check 
	for (int i = 0;i < n;i++) 
	    if      (REF.length() == 1 && alts[i].length() == 1) ; // SNP
	    else if (REF.length() == 0 || alts[i].length() == 0) {
		System.err.println("Skipping invalid record:");
		System.err.println(" " + line);
		return null;
	    }
	
	// Adjustment of first base
	if (REF.length() > 0) {
	    boolean same = true;
	    for (int i = 0;i < n;i++) 
		if (alts[i].length() == 0 ||
		    REF.charAt(0) != alts[i].charAt(0)) {
		    same = false;
		    break;
		}
	    if (same) {
		pos++;
		REF = REF.substring(1);
		for (int i = 0;i < n;i++) alts[i] = alts[i].substring(1);
	    }
	}

	// Adjustment of last
	if (REF.length() > 0) {
	    boolean same = true;
	    int indREF = REF.length() - 1;
	    for (int i = 0;i < n;i++) {
		int len = alts[i].length();
		if (len == 0 || 
		    REF.charAt(indREF) != alts[i].charAt(len - 1)) {
		    same = false;
		    break;
		}
	    }
	    if (same) {
		REF = REF.substring(0,indREF);
		for (int i = 0;i < n;i++)
		    alts[i] = alts[i].substring(0,alts[i].length() - 1);
	    }

	    // System.out.print(REF);
	    // for (int i = 0;i < n;i++)
	    //     System.out.print(" " + alts[i]);
	    // System.out.println(pos);
	}

	return new Variant(chr,pos,REF.length(),alts,phase);
    }

    private void readLine()
    {
	_line = null;
	try {
	    _line = _br.readLine();
	} catch (Exception ex) { }

	if (_line == null) {
	    try {
		_br.close();
		_fin.close();
	    } catch (Exception ex) { }
	    _br  = null;
	    _fin = null;
	}
    }

    private int getGenotypeIndex(String record)
    {
	StringTokenizer words = new StringTokenizer(record,":");
	int ret = 0;
	while (words.hasMoreTokens())
	    if (words.nextToken().equals("GT")) return ret;
	    else                                ret++;
	return -1;
    }

    public static int getChromIndex(String chr)
    {
	int ret = -1;
	String tmp = chr;
	if (chr.length() > 3 &&
	    chr.substring(0,3).equalsIgnoreCase("chr"))
	    chr = chr.substring(3);
	if      (chr.equalsIgnoreCase("X"))  ret = 23;
	else if (chr.equalsIgnoreCase("Y"))  ret = 24;
	else if (chr.equalsIgnoreCase("M") ||
		 chr.equalsIgnoreCase("MT")) ret = 25;
	// Note: Contigs added for hs37d5 compatibility.
    	else if (chr.equalsIgnoreCase("GL000207.1")) ret = 26;
    	else if (chr.equalsIgnoreCase("GL000226.1")) ret = 27;
    	else if (chr.equalsIgnoreCase("GL000229.1")) ret = 28;
    	else if (chr.equalsIgnoreCase("GL000231.1")) ret = 29;
    	else if (chr.equalsIgnoreCase("GL000210.1")) ret = 30;
    	else if (chr.equalsIgnoreCase("GL000239.1")) ret = 31;
    	else if (chr.equalsIgnoreCase("GL000235.1")) ret = 32;
    	else if (chr.equalsIgnoreCase("GL000201.1")) ret = 33;
    	else if (chr.equalsIgnoreCase("GL000247.1")) ret = 34;
    	else if (chr.equalsIgnoreCase("GL000245.1")) ret = 35;
    	else if (chr.equalsIgnoreCase("GL000197.1")) ret = 36;
    	else if (chr.equalsIgnoreCase("GL000203.1")) ret = 37;
    	else if (chr.equalsIgnoreCase("GL000246.1")) ret = 38;
    	else if (chr.equalsIgnoreCase("GL000249.1")) ret = 39;
    	else if (chr.equalsIgnoreCase("GL000196.1")) ret = 40;
    	else if (chr.equalsIgnoreCase("GL000248.1")) ret = 41;
    	else if (chr.equalsIgnoreCase("GL000244.1")) ret = 42;
    	else if (chr.equalsIgnoreCase("GL000238.1")) ret = 43;
    	else if (chr.equalsIgnoreCase("GL000202.1")) ret = 44;
    	else if (chr.equalsIgnoreCase("GL000234.1")) ret = 45;
    	else if (chr.equalsIgnoreCase("GL000232.1")) ret = 46;
    	else if (chr.equalsIgnoreCase("GL000206.1")) ret = 47;
    	else if (chr.equalsIgnoreCase("GL000240.1")) ret = 48;
    	else if (chr.equalsIgnoreCase("GL000236.1")) ret = 49;
    	else if (chr.equalsIgnoreCase("GL000241.1")) ret = 50;
    	else if (chr.equalsIgnoreCase("GL000243.1")) ret = 51;
    	else if (chr.equalsIgnoreCase("GL000242.1")) ret = 52;
    	else if (chr.equalsIgnoreCase("GL000230.1")) ret = 53;
    	else if (chr.equalsIgnoreCase("GL000237.1")) ret = 54;
    	else if (chr.equalsIgnoreCase("GL000233.1")) ret = 55;
    	else if (chr.equalsIgnoreCase("GL000204.1")) ret = 56;
    	else if (chr.equalsIgnoreCase("GL000198.1")) ret = 57;
    	else if (chr.equalsIgnoreCase("GL000208.1")) ret = 58;
    	else if (chr.equalsIgnoreCase("GL000191.1")) ret = 59;
    	else if (chr.equalsIgnoreCase("GL000227.1")) ret = 60;
    	else if (chr.equalsIgnoreCase("GL000228.1")) ret = 61;
    	else if (chr.equalsIgnoreCase("GL000214.1")) ret = 62;
    	else if (chr.equalsIgnoreCase("GL000221.1")) ret = 63;
    	else if (chr.equalsIgnoreCase("GL000209.1")) ret = 64;
    	else if (chr.equalsIgnoreCase("GL000218.1")) ret = 65;
    	else if (chr.equalsIgnoreCase("GL000220.1")) ret = 66;
    	else if (chr.equalsIgnoreCase("GL000213.1")) ret = 67;
    	else if (chr.equalsIgnoreCase("GL000211.1")) ret = 68;
    	else if (chr.equalsIgnoreCase("GL000199.1")) ret = 69;
    	else if (chr.equalsIgnoreCase("GL000217.1")) ret = 70;
    	else if (chr.equalsIgnoreCase("GL000216.1")) ret = 71;
    	else if (chr.equalsIgnoreCase("GL000215.1")) ret = 72;
    	else if (chr.equalsIgnoreCase("GL000205.1")) ret = 73;
    	else if (chr.equalsIgnoreCase("GL000219.1")) ret = 74;
    	else if (chr.equalsIgnoreCase("GL000224.1")) ret = 75;
    	else if (chr.equalsIgnoreCase("GL000223.1")) ret = 76;
    	else if (chr.equalsIgnoreCase("GL000195.1")) ret = 77;
    	else if (chr.equalsIgnoreCase("GL000212.1")) ret = 78;
    	else if (chr.equalsIgnoreCase("GL000222.1")) ret = 79;
    	else if (chr.equalsIgnoreCase("GL000200.1")) ret = 80;
    	else if (chr.equalsIgnoreCase("GL000193.1")) ret = 81;
    	else if (chr.equalsIgnoreCase("GL000194.1")) ret = 82;
    	else if (chr.equalsIgnoreCase("GL000225.1")) ret = 83;
    	else if (chr.equalsIgnoreCase("GL000192.1")) ret = 84;
    	else if (chr.equalsIgnoreCase("NC_007605")) ret = 85;
    	else if (chr.equalsIgnoreCase("hs37d5")) ret = 86;
	else
	    try {
		ret = Integer.parseInt(chr);
	    } catch (Exception e) {
		System.err.println("Unknown chromosome " + chr + ".");
	    }
	
	return ret;
    }
}
