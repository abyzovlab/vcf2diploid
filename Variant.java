//--- Java imports ---
import java.util.*;
import java.io.*;
import java.lang.reflect.*;

class Variant
{
    private int      _chr = -1, _pos = -1, _del = -1;
    private String[] _alts;
    private int      _maternal = 0, _paternal = 0;
    private boolean  _isPhased = false; // Phasing
    private static final Random _rand = new Random();
    
    public Variant(int chr,int pos,int del,String[] alts,String phase)
    {
	_chr  = chr;
	_pos  = pos;
	_del  = del;
	_alts = new String[alts.length];
	for (int i = 0;i < alts.length;i++) _alts[i] = new String(alts[i]);
	
	phase = phase.trim();
	boolean strangePhase = false;
	if (phase.length() == 1) {
	    if (Character.isDigit(phase.charAt(0))) {
		int val = Integer.parseInt(phase);
		if (chr == 22) {
		    _maternal = val;
		    _isPhased = true;
		} else if (chr == 23) {
		    _paternal = val;
		    _isPhased = true;
		} else strangePhase = true;
	    } else strangePhase = true;
	} else if (phase.length() == 3) {
	    char c1      = phase.charAt(0);
	    char c2      = phase.charAt(2);
	    char phasing = phase.charAt(1);

	    if (Character.isDigit(c1) && Character.isDigit(c2)) {
		_paternal = Character.digit(c1,10);
		_maternal = Character.digit(c2,10);
		if (phasing == '|') _isPhased = true;
	    } else strangePhase = true;
	} else strangePhase = true;

	if (strangePhase)
	    System.err.println("Unreconized phasing '" + phase + "'.");
    }
    
    public int    chromosome() { return _chr; }
    public int    position()   { return _pos; }
    public int    deletion()   { return _del; }
    public int    maternal()   { return _maternal; }
    public int    paternal()   { return _paternal; }
    public String insertion(int ind)
    {
	if (ind <= 0 || ind > _alts.length) return "";
	return _alts[ind - 1];
    }
    public int    variantBases()
    {
	int ret = _del;
	for (int i = 0;i < _alts.length;i++)
	    if (_del != _alts[i].length()) ret += _alts[i].length();
	return ret;
    }
    public boolean isPhased() { return _isPhased; }
    public void randomizeHaplotype()
    {
	if (_rand.nextDouble() > 0.5) return;
	int tmp   = _paternal;
	_paternal = _maternal;
	_maternal = tmp;
	return;
    }

}
