//--- Java imports ---
import java.util.*;
import java.io.*;
import java.lang.reflect.*;

/**
 * Class to construct diploid genome from genome reference and genome variants
 * calls in VCF format.
 *
 * @author Alexej Abyzov
 */
public class VCF2diploid
{
    private final static char DELETED_BASE = '~';

    private String[] _chrFiles = null, _vcfFiles = null;
    private String   _id = "";
    private ArrayList<Variant>[] _variants = new ArrayList[25];
	
    public VCF2diploid(String[] chrFiles,String[] vcfFiles,
		       String id,boolean pass)
    {
	_chrFiles = chrFiles;
	_vcfFiles = vcfFiles;
	if (id != null) _id = id;

	for (int i = 0;i < _variants.length;i++)
	    _variants[i] = new ArrayList<Variant>(128);

	for (int i = 0;i < vcfFiles.length;i++) {
	    VCFparser parser = new VCFparser(vcfFiles[i],_id,pass);
	    int n_ev = 0,var_nucs = 0;
	    while (parser.hasMoreInput()) {
	    	Variant var = parser.parseLine();
	    	if (var == null) continue;
		if (var.maternal() == 0 &&
		    var.paternal() == 0) continue;
	    	int chr = var.chromosome();
	    	if (chr <= 0 || chr > _variants.length) continue;
	    	_variants[chr - 1].add(var);
		n_ev++;
		var_nucs += var.variantBases();
	    }
	    System.out.println(vcfFiles[i] + ": " + n_ev + " variants, " +
			       var_nucs + " variant bases");
	}

	// for (int i = 0;i < _variants.length;i++)
	//     System.out.println((i + 1) + " " + _variants[i].size());
    }

    public void makeDiploid()
    {
	StringBuffer paternal_chains = new StringBuffer("");
	StringBuffer maternal_chains = new StringBuffer("");
	int chain_id = 1;
	for (int f = 0;f < _chrFiles.length;f++) {
	    // System.gc();
	    // System.out.println("Before parsing sequence.");
	    // System.out.println("Used memory = " +
	    // 		       (Runtime.getRuntime().totalMemory() -
	    // 			Runtime.getRuntime().freeMemory()));
	    Sequence[] seqs = parseSequences(_chrFiles[f]);
	    // System.gc();
	    // System.out.println("After parsing sequence.");
	    // System.out.println("Used memory = " +
	    // 		       (Runtime.getRuntime().totalMemory() -
	    // 			Runtime.getRuntime().freeMemory()));
	    for (int s = 0;s < seqs.length;s++) {
		System.out.println("Working on " + seqs[s].getName() + "...");
		int index = VCFparser.getChromIndex(seqs[s].getName());
		if (index <= 0 || index > _variants.length) continue;
		ArrayList<Variant> varList = _variants[index - 1];

		Sequence ref_seq = seqs[s];

		int len = ref_seq.length();
		byte[] maternal_seq = new byte[len];
		byte[] paternal_seq = new byte[len];
		byte[] ins_flag     = new byte[len];
		// Flag specification:
		// b -- insertion in both haplotypes
		// p -- insertion in paternal haplotype
		// m -- insertion in maternal haplotype

		for (int c = 0;c < len;c++)
		    maternal_seq[c] = paternal_seq[c] = ref_seq.byteAt(c);

		Hashtable<Integer,String> pat_ins_seq =
		    new Hashtable<Integer,String>(150);
		Hashtable<Integer,String> mat_ins_seq =
		    new Hashtable<Integer,String>(150);

		int n_var_pat  = 0,n_var_mat  = 0;
		int n_base_pat = 0,n_base_mat = 0;
		ListIterator<Variant> it = varList.listIterator();
		while (it.hasNext()) {
		    Variant var = it.next();
		    int     pos = var.position() - 1;
		    int     del = var.deletion();
		    if (!var.isPhased()) var.randomizeHaplotype();
		    if (var.paternal() > 0)
			if (addVariant(paternal_seq,ref_seq,
				       pos,del,var.insertion(var.paternal()),
				       pat_ins_seq)) {
			    n_var_pat++;
			    n_base_pat += var.variantBases();
			}
		    if (var.maternal() > 0)
			if (addVariant(maternal_seq,ref_seq,
				       pos,del,var.insertion(var.maternal()),
				       mat_ins_seq)) {
			    n_var_mat++;
			    n_base_mat += var.variantBases();
			}
		}

		writeMap(seqs[s],paternal_seq,maternal_seq,
			 pat_ins_seq,mat_ins_seq);
		writeDiploid(ref_seq,paternal_seq,maternal_seq,
			     pat_ins_seq,mat_ins_seq);

		paternal_chains.append(makeChains(seqs[s].getName(),
						  paternalName(seqs[s].getName()),
						  paternal_seq,pat_ins_seq,
						  chain_id));
		maternal_chains.append(makeChains(seqs[s].getName(),
						  maternalName(seqs[s].getName()),
						  maternal_seq,mat_ins_seq,
						  chain_id));
		chain_id++;
		
		System.out.println("Applied " + n_var_pat + " variants " +
				   n_base_pat + " bases to " +
				   "paternal genome.");
		System.out.println("Applied " + n_var_mat + " variants " +
				   n_base_mat + " bases to " +
				   "maternal genome.");
	    }
	}

	try {
	    FileWriter fw = new FileWriter(new File("paternal.chain"));
	    BufferedWriter bw = new BufferedWriter(fw);
	    bw.write(paternal_chains.toString());
	    bw.newLine();
	    bw.close();
	    fw.close();
	} catch (Exception ex) {
	    System.err.println(ex.toString());
	}

	try {
	    FileWriter fw = new FileWriter(new File("maternal.chain"));
	    BufferedWriter bw = new BufferedWriter(fw);
	    bw.write(maternal_chains.toString());
	    bw.newLine();
	    bw.close();
	    fw.close();
	} catch (Exception ex) {
	    System.err.println(ex.toString());
	}
    }

    private boolean addVariant(byte[] new_seq,Sequence ref_seq,
			       int pos,int del,String ins,
			       Hashtable<Integer,String> ins_seq)
    {
	boolean overlap = false;

	if (pos >= new_seq.length || pos + del >= new_seq.length) {
	    System.err.println("Variant out of chromosome bounds at " +
			       ref_seq.getName() + ":" + (pos + 1) +
			       ", (del,ins) of (" +
			       del + "," + ins + ").");
	    System.err.println("Skipping.");
	    return false;
	}

	for (int p = pos;p < pos + del;p++)
	    if (new_seq[p] == DELETED_BASE ||
		new_seq[p] != ref_seq.byteAt(p)) overlap = true;
	
	if (overlap) {
	    System.err.println("Variant overlap at " +
			       ref_seq.getName() + ":" + (pos + 1) +
			       ", (del,ins) of (" +
			       del + "," + ins + ").");
	    System.err.println("Skipping.");
	    return false;
	}

	if (del == 1 && ins.length() == 1) { // SNP
	    if (Character.isLowerCase((char)ref_seq.byteAt(pos)))
		new_seq[pos] = (byte)Character.toLowerCase(ins.charAt(0));
	    else
		new_seq[pos] = (byte)Character.toUpperCase(ins.charAt(0));
	} else { // Indel, SV
	    //if (del <= 0) pos++; // Position adjustment
	    if (ins_seq.get(pos) != null) {
		System.err.println("Multiple insertions at " + 
				   ref_seq.getName() + ":" + (pos + 1));
		System.err.println("Skipping variant with (del,ins) of (" +
				   del + "," + ins + ").");
		return false;
	    }
	    for (int p = pos;p < pos + del;p++) 
		new_seq[p] = DELETED_BASE;

	    // if (for_pat && for_mat) ins_flag[pos] = 'b';
	    // else if (for_pat)       ins_flag[pos] = 'p';
	    // else if (for_mat)       ins_flag[pos] = 'm';

	    if (ins.length() > 0) ins_seq.put(new Integer(pos),ins);
	}

	return true;
    }

    private String makeChains(String ref_name,String der_name,
			      byte[] genome,Hashtable<Integer,String> ins_seq,
			      int id)
    {
	boolean[] ins_flag = new boolean[genome.length];
	Enumeration enm = ins_seq.keys();
	while (enm.hasMoreElements()) {
	    Integer key = (Integer)enm.nextElement();
	    ins_flag[key.intValue()] = true;
	}
	int ref_len = genome.length;
	int der_len = 0;
	int score   = 0;
	for (int p = 0;p < genome.length;p++) {
	    if (ins_flag[p]) der_len += ins_seq.get(p).length();
	    if (genome[p] != DELETED_BASE) {
		der_len++;
		score++;
	    }
	}

	StringWriter ret = new StringWriter();
        PrintWriter  wr  = new PrintWriter(ret);
	wr.println("chain " + score + " " +
		   ref_name + " " + ref_len + " + 0 " + ref_len + " " +
		   der_name + " " + der_len + " + 0 " + der_len + " " + id);
	int size = 0, dref = 0, dder = 0;
	boolean flag = false;
	for (int p = 0;p < genome.length;p++) {
	    if (ins_flag[p]) {
		dder += ins_seq.get(p).length();
		flag = true;
	    }
	    if (genome[p] == DELETED_BASE) {
		dref++;
		flag = true;
	    } else { // Normal base
		if (flag) {
		    wr.println(size + " " + dref + " " + dder);
		    size = dref = dder = 0;
		    flag = false;
		}
		size++;
	    }
	}
	wr.println(size);
	wr.println();

	return ret.toString();
    }

    private void writeMap(Sequence ref_seq,
			  byte[] paternal,byte[] maternal,
			  Hashtable<Integer,String> ins_seq_pat,
			  Hashtable<Integer,String> ins_seq_mat)
    {
	if (paternal.length != maternal.length) {
	    System.err.println("Paternal and maternal genomes are of " +
			       "different lengths. Making map aborted.");
	    return;
	}

	boolean[] ins_flag_pat = new boolean[paternal.length];
	Enumeration enm = ins_seq_pat.keys();
	while (enm.hasMoreElements()) {
	    Integer key = (Integer)enm.nextElement();
	    ins_flag_pat[key.intValue()] = true;
	}
	boolean[] ins_flag_mat = new boolean[maternal.length];
	enm = ins_seq_mat.keys();
	while (enm.hasMoreElements()) {
	    Integer key = (Integer)enm.nextElement();
	    ins_flag_mat[key.intValue()] = true;
	}


	int NOT_IN_GENOME = 0;
	String file_name = ref_seq.getName() + "_" + _id + ".map";
	try {
	    FileWriter fw = new FileWriter(new File(file_name));
	    BufferedWriter bw = new BufferedWriter(fw);
	    bw.write("#REF\tPAT\tMAT");
	    bw.newLine();
	    int ri =             1, pi =             1, mi =             1;
	    int pr = NOT_IN_GENOME, pp = NOT_IN_GENOME, pm = NOT_IN_GENOME;
	    for (int p = 0;p < paternal.length;p++) {
		boolean for_pat  = ins_flag_pat[p];
		boolean for_mat  = ins_flag_mat[p];
		boolean for_both = for_pat && for_mat;
		if (for_both)
		    for_both = ins_seq_pat.get(p).equals(ins_seq_mat.get(p));
		if (for_pat || for_mat || for_both)
		    if (pr > 0 || pp > 0 || pm > 0) {
			bw.write(pr + "\t" + pp + "\t" + pm);
			bw.newLine();
			pr = pp = pm = 0;
		    }
		if (for_both) {
		    String new_seq = ins_seq_pat.get(p);
		    if (new_seq != null) {
			bw.write(NOT_IN_GENOME + "\t" + pi + "\t" + mi);
			bw.newLine();
			for (int i = 0;i < new_seq.length();i++) {
			    pi++; mi++;
			}
		    }
		} else {
		    if (for_pat) {
			String new_seq = ins_seq_pat.get(p);
			if (new_seq != null) {
			    bw.write(NOT_IN_GENOME + "\t" + pi + "\t" +
				     NOT_IN_GENOME);
			    bw.newLine();
			    for (int i = 0;i < new_seq.length();i++) pi++;
			}
		    }
		    if (for_mat) {
			String new_seq = ins_seq_mat.get(p);
			if (new_seq != null) {
			    bw.write(NOT_IN_GENOME + "\t" + NOT_IN_GENOME +
				     "\t" + mi);
			    bw.newLine();
			    for (int i = 0;i < new_seq.length();i++) mi++;
			}
		    }
		}
		int ref = ri++,pat = NOT_IN_GENOME,mat = NOT_IN_GENOME;
		if (paternal[p] != DELETED_BASE) pat = pi++;
		if (maternal[p] != DELETED_BASE) mat = mi++;
		if (pr == NOT_IN_GENOME &&
		    pp == NOT_IN_GENOME &&
		    pm == NOT_IN_GENOME) { // Initiation
		    pr = ref; pp = pat; pm = mat;
		} else {
		    boolean cand = 
			((pat == NOT_IN_GENOME && pp == NOT_IN_GENOME) ||
			 ref - pr == pat - pp) &&
			((mat == NOT_IN_GENOME && pm == NOT_IN_GENOME) ||
			 ref - pr == mat - pm);
		    if (!cand) {
			bw.write(pr + "\t" + pp + "\t" + pm);
			bw.newLine();
			pr = ref; pp = pat; pm = mat;
		    }
		}
	    }
	    if (pr > 0 || pp > 0 || pm > 0) {
		bw.write(pr + "\t" + pp + "\t" + pm);
		bw.newLine();
	    }
	    bw.close();
	    fw.close();
	} catch (Exception ex) {
	    System.err.println(ex.toString());
	}
    }

    private void writeDiploid(Sequence ref_seq,
			      byte[] paternal,byte[] maternal,
			      Hashtable<Integer,String> pat_ins_seq,
			      Hashtable<Integer,String> mat_ins_seq)
    {
	String file_name = paternalName(ref_seq.getName() + "_" + _id) + ".fa";
	String name      = paternalName(ref_seq.getName());
	try {
	    FileWriter fw = new FileWriter(new File(file_name));
	    BufferedWriter bw = new BufferedWriter(fw);
	    writeGenome(bw,name,paternal,pat_ins_seq);
	    bw.close();
	    fw.close();
	} catch (Exception ex) {
	    System.err.println(ex.toString());
	}

	file_name = maternalName(ref_seq.getName() + "_" + _id) + ".fa";
	name      = maternalName(ref_seq.getName());
	try {
	    FileWriter fw = new FileWriter(new File(file_name));
	    BufferedWriter bw = new BufferedWriter(fw);
	    writeGenome(bw,name,maternal,mat_ins_seq);
	    bw.close();
	    fw.close();
	} catch (Exception ex) {
	    System.err.println(ex.toString());
	}
    }

    private void writeGenome(BufferedWriter bw,String name,
			     byte[] genome,Hashtable<Integer,String> ins_seq)
	throws Exception
    {
	final int line_width = 50;

	boolean[] ins_flag = new boolean[genome.length];
	Enumeration enm = ins_seq.keys();
	while (enm.hasMoreElements()) {
	    Integer key = (Integer)enm.nextElement();
	    ins_flag[key.intValue()] = true;
	}

	bw.write(">" + name);
	bw.newLine();

	StringBuffer line = new StringBuffer();
	for (int p = 0;p < genome.length;p++) {
	    if (ins_flag[p]) {
		String new_seq = ins_seq.get(p);
		if (new_seq != null && new_seq.length() > 0)
		    line.append(new_seq);
	    }
	    if (genome[p] != DELETED_BASE) line.append((char)genome[p]);
	    while (line.length() >= line_width) {
		bw.write(line.toString(),0,line_width);
		bw.newLine();
		line.delete(0,line_width);
	    }
	}
	while (line.length() > 0) {
	    int n = line.length();
	    if (line_width < n) n = line_width; 
	    bw.write(line.toString(),0,n);
	    bw.newLine();
	    line.delete(0,n);
	}
    }

    private Sequence[] parseSequences(String fileName)
    {
	ArrayList<Sequence> ss = new ArrayList<Sequence>(1);
	String header = "";
	File file = new File(fileName);
	byte[] seq = new byte[300000000];
	int index = 0;
	try {
	    FileReader     fr = new FileReader(file);
	    BufferedReader br = new BufferedReader(fr);
	    String line = new String("");
	    while ((line = br.readLine()) != null) {
		if (line.startsWith(">")) {
		    if (index > 0)
			ss.add(new Sequence(header,seq,index));
		    header = line;
		    index  = 0;
		} else {
		    for (int i = 0;i < line.length();i++)
			seq[index++] = (byte)line.charAt(i);
		}
	    }
	    br.close();
	    fr.close();
	} catch (Exception ex) {
	    System.err.println("Can't open file " + fileName);
	    System.err.println(ex.toString());
	}
	if (index > 0) ss.add(new Sequence(header,seq,index));
	return ss.toArray(new Sequence[0]);
    }

    private String maternalName(String name)
    {
	return (name + "_maternal");
    }

    private String paternalName(String name)
    {
	return (name + "_paternal");
    }
    
    /**
     * Main function.
     */
    public static void main(String[] args)
    {
	String VERSION = "vcf2diploid - v0.2.6";

	ArrayList<String> chrFiles = new ArrayList<String>(1);
	ArrayList<String> vcfFiles = new ArrayList<String>(1);
	String id = "";
	boolean pass = false;

	String usage = "Usage:\n";
	usage += "\tvcf2diploid -id sample_id [-pass] ";
	usage += "-chr file.fa ... ";
	usage += "[-vcf file.vcf ...]\n";
	usage += "\tvcf2diploid -version\n";

	for (int i = 0;i < args.length;i++) {
	    if (args[i].equals("-vcf")) {
		while (++i < args.length)
		    if (args[i].charAt(0) != '-') vcfFiles.add(args[i]);
		    else { i--;	break; }
	    } else if (args[i].equals("-chr")) {
		while (++i < args.length) 
		    if (args[i].charAt(0) != '-') chrFiles.add(args[i]);
		    else { i--; break; }
	    } else if (args[i].equals("-id")) {
		if (++i < args.length) id = args[i];
	    } else if (args[i].equals("-version")) {
		System.out.println(VERSION);
		return;
	    } else if (args[i].equals("-pass")) {
		pass = true;
	    }
	}

	if (id.length() <= 0) {
	    System.err.println("No sample id is given.\n");
	    System.err.println(usage);
	    return;
	}

	if (chrFiles.size() == 0) {
	    System.err.println("No chromosome file(s) is given!\n");
	    System.err.println(usage);
	    return;
	}

	if (vcfFiles.size() == 0)
	    System.err.println("No VCF file(s) is given!");

	VCF2diploid maker =
	    new VCF2diploid(chrFiles.toArray(new String[0]),
			    vcfFiles.toArray(new String[0]),
			    id,pass);
	maker.makeDiploid();
    }


    private class Sequence
    {
	private String _header = "",_name = "";
	private byte[] _seq = null;

	public Sequence(String header,byte[] seq,int len)
	{
	    _header = header;
	    if (len < 0) len = 0;
	    if (seq.length < len) len = seq.length;
	    _seq = new byte[len];
	    for (int i = 0;i < len;i++) _seq[i] = seq[i];
	    StringTokenizer toks = new StringTokenizer(header.substring(1));
	    if (toks.hasMoreTokens()) _name = toks.nextToken();
	}

	public String getName()     { return _name; }
	public String getHeader()   { return _header; }
	public int    length()      { return _seq.length; }
	public byte   byteAt(int p) { return _seq[p]; }
    }
}