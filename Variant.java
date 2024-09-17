//--- Java imports ---
import java.util.*;
import java.io.*;
import java.lang.reflect.*;

class Variant {
	private String chr;
	private int _pos = -1;
	private int _del = -1;
	private String[] _alts;
	private int _maternal = 0, _paternal = 0;
	private boolean _isPhased = false; // Phasing
	private static final Random _rand = new Random();

	public Variant(String chr,int pos,int del,String[] alts,String phase) {
		this.chr  = chr;
		_pos  = pos;
		_del  = del;
		_alts = new String[alts.length];
		for (int i = 0;i < alts.length;i++) _alts[i] = new String(alts[i]);

		phase = phase.trim();
		boolean strangePhase = false;
		if (phase.length() == 1) {
			if (Character.isDigit(phase.charAt(0))) {
				int val = Integer.parseInt(phase);
				if (chr.contains("chrX")) {
					_maternal = val;
					_isPhased = true;
				} else if (chr.contains("chrY")) {
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

	public String chromosome() { return chr; }
	public int position() { return _pos; }
	public int deletion() { return _del; }
	public int maternal() { return _maternal; }
	public int paternal() { return _paternal; }
	public String insertion(int ind) {
		if (ind <= 0 || ind > _alts.length) return "";
		return _alts[ind - 1];
	}
	public int variantBases() {
		int ret = _del;
		for (int i = 0;i < _alts.length;i++)
			if (_del != _alts[i].length()) ret += _alts[i].length();
		return ret;
	}
	public boolean isPhased() { return _isPhased; }
	public void randomizeHaplotype() {
		if (_rand.nextDouble() > 0.5) return;
		int tmp   = _paternal;
		_paternal = _maternal;
		_maternal = tmp;
		return;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(chr);
		sb.append("\t");
		sb.append(_pos);
		sb.append("\t");
		sb.append(_del);
		sb.append("\t");
		sb.append(_alts.length);
		sb.append("\t");
		sb.append(_maternal);
		sb.append("\t");
		sb.append(_paternal);
		sb.append("\t");
		sb.append(_isPhased);
		return sb.toString();
	}
}
