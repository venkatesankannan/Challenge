import java.util.Arrays;
import java.util.ArrayList;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayDeque;
import java.io.BufferedReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.Collections;
import java.io.FileReader;
import java.io.IOException;
import java.util.AbstractCollection;
import java.util.StringTokenizer;
import java.util.Set;


public class Solution {
	public static void main(String[] args) {
		InputStream inputStream = System.in;
		OutputStream outputStream = System.out;
		InputReader in = new InputReader(inputStream);
		OutputWriter out = new OutputWriter(outputStream);
		CancerIQ solver = new CancerIQ();
		solver.solve(1, in, out);
		out.close();
	}
}

class CancerIQ {

	static class Range {
		final int N;
		final int root;
		int[] rangeSum;
		int[] sum;
		int[] rangeMax;

		public Range(int N) {
			this.N = N;
			root = (int)Math.sqrt(N);
			rangeSum = new int[getIndex(N - 1) + 1];
			rangeMax = new int[rangeSum.length];
			sum = new int[N];
		}

		int getSum(int queryStart, int queryEnd) {
			int ans = 0;

			int first = getIndex(queryStart);
			if (queryStart != getStart(queryStart)) {
				first++;
			}

			int last = getIndex(queryEnd);
			if (queryEnd != getEnd(queryEnd)) {
				last--;
			}

			if (first <= last) {
				for (int i = queryStart; i < first * root; i++) {
					ans += sum[i];
				}

				for (int i = first; i <= last; i++) {
					ans += rangeSum[i];
				}

				for (int i = root * (last + 1); i <= queryEnd; i++) {
					ans += sum[i];
				}
			} else {
				for (int i = queryStart; i <= queryEnd; i++) {
					ans += sum[i];
				}
			}


			return ans;
		}

		int getMax(int queryStart, int queryEnd) {
			int first = getIndex(queryStart);
			if (queryStart != getStart(queryStart)) {
				first++;
			}

			int last = getIndex(queryEnd);
			if (queryEnd != getEnd(queryEnd)) {
				last--;
			}

			int max = 0;

			if (first <= last) {
				for (int i = queryEnd; i >= root * (last + 1); i--) {
					max = Math.max(max + sum[i], sum[i]);
				}

				for (int i = last; i >= first; i--) {
					max = Math.max(max + rangeSum[i], rangeMax[i]);
				}

				for (int i = first * root - 1; i >= queryStart; i--) {
					max = Math.max(max + sum[i], sum[i]);
				}
			} else {
				for (int i = queryEnd; i >= queryStart; i--) {
					max = Math.max(max + sum[i], sum[i]);
				}
			}

			return max;
		}

		void add(int pos, int value) {
			int range = getIndex(pos);

			sum[pos] += value;

			rangeSum[range] += value;

			int start = getStart(pos);
			int end = getEnd(pos);

			rangeMax[range] = 0;

			for (int i = end; i >= start; i--) {
				rangeMax[range] = Math.max(rangeMax[range] + sum[i], sum[i]);
			}
		}

		int getIndex(int pos) {
			return pos / root;
		}

		int getStart(int pos) {
			return root * getIndex(pos);
		}

		int getEnd(int pos) {
			return Math.min(getStart(pos) + root, N) - 1;
		}
	}

	int N;
	int[] parent;
	int[] chain;
	int[] index;
	ArrayList<ArrayList<Integer>> chains;
	ArrayList<ArrayList<Integer>> graph = new ArrayList<>();
	Range[] ranges;

	public void solve(int testNumber, InputReader in, OutputWriter out) {
		N = in.readInt();

		for (int i = 0; i < N; i++) {
			graph.add(new ArrayList<>());
		}

		for (int i = 1; i < N; i++) {
			int u = in.readInt() - 1;
			int v = in.readInt() - 1;
			graph.get(u).add(v);
			graph.get(v).add(u);
		}

		preprocess();

		int Q = in.readInt();

		while (Q-- > 0) {
			if (in.next().equals("add")) {
				int u = in.readInt() - 1;
				int value = in.readInt();
				add(u, value);
				
			} else {
				int a = in.readInt() - 1;
				int b = in.readInt() - 1;
				int lca = getLca(a, b);
				int maxA = getMax(lca, a, 0);
				int maxB = getMax(lca, b, 0);
				int max = Math.max(maxA, maxB);
				if (lca != 0) {
					
					max += getSum(parent[lca]);
				}
				out.println(max);
			}
			
		}
	}

	void add(int node, int value) {
		ranges[chain[node]].add(index[node], value);
	}

	int getMax(int lca, int a, int max) {
		int start = chain[lca] == chain[a] ? index[lca] : 0;
		int end = index[a];
		int sum = ranges[chain[a]].getSum(start, end);
		max = Math.max(max + sum, ranges[chain[a]].getMax(start, end));
		if (chain[lca] == chain[a]) {
			return max;
		}
		return getMax(lca, parent[chains.get(chain[a]).get(0)], max);
	}

	int getSum(int a) {
		int sum = ranges[chain[a]].getSum(0, index[a]);
		if (chain[a] == 0) {
			return sum;
		}
		return sum + getSum(parent[chains.get(chain[a]).get(0)]);
	}

	int getLca(int a, int b) {
		if (chain[a] == chain[b]) {
			return chains.get(chain[a]).get(Math.min(index[a], index[b]));
		}
		if (chain[a] < chain[b]) {
			a ^= b;
			b ^= a;
			a ^= b;
		}
		return getLca(parent[chains.get(chain[a]).get(0)], b);
	}

	void preprocess() {
		Tree.HLD tempInfo = Tree.getHLD(graph);
		parent = tempInfo.parent;
		chain = tempInfo.chain;
		chains = tempInfo.chains;
		index = tempInfo.index;
		ranges = new Range[chains.size()];
		for (int i = 0; i < chains.size(); i++) {
			ranges[i] = new Range(chains.get(i).size());
		}
	}

}

class InputReader {
	private BufferedReader input;
	private StringTokenizer line = new StringTokenizer("");

	public InputReader(InputStream in) {
		input = new BufferedReader(new InputStreamReader(in));
	}

	public void fill() {
		try {
			if(!line.hasMoreTokens()) line = new StringTokenizer(input.readLine());
		} catch(IOException io) { io.printStackTrace(); System.exit(0);}
	}

	public String next() {
		fill();
		return line.nextToken();
	}

	public int readInt() {
		fill();
		return Integer.parseInt(line.nextToken());
	}

}

class OutputWriter {
	private PrintWriter output;

	public OutputWriter(OutputStream out) {
		output = new PrintWriter(out);
	}

	public void println(Object o) {
		output.println(o);
	}

	public void close() {
		output.close();
	}
}

class Tree {

	public static class HLD {
		public int[] parent;
		public int[] chain;
		public int[] index;

		public ArrayList<ArrayList<Integer>> chains;

		private HLD(int N) {
			parent = new int[N];
			chain = new int[N];
			index = new int[N];
			chains = new ArrayList<>();
		}

		public String toString() {
			String ans = chains.toString();
			for (int i = 0; i < parent.length; i++) {
				ans += System.lineSeparator();
				ans += String.format("node: %d, parent: %d, chain: %d, index: %d",
						i, parent[i], chain[i], index[i]);
			}
			return ans;
		}
	}

	public static HLD getHLD(ArrayList<ArrayList<Integer>> graph) {
		return getHLD(graph, 0);
	}

	public static HLD getHLD(ArrayList<ArrayList<Integer>> graph, int root) {
		final int N = graph.size();

		HLD info = new HLD(N);

		Arrays.fill(info.parent, -1);

		ArrayList<Integer> order = new ArrayList<>();

		ArrayDeque<Integer> stack = new ArrayDeque<>();
		stack.push(root);
		while (!stack.isEmpty()) {
			int node = stack.pop();
			order.add(node);
			for (int child: graph.get(node)) if (child != info.parent[node]) {
				info.parent[child] = node;
				stack.add(child);
			}
		}

		int[] size = new int[N];

		for (int i = N - 1; i >= 0; i--) {
			int node = order.get(i);
			size[node] = 1;
			for (int child: graph.get(node)) if (child != info.parent[node]) {
				size[node] += size[child];
			}
		}

		info.chains.add(new ArrayList<>(Collections.singleton(root)));
		info.chain[root] = 0;
		info.index[root] = 0;

		for (int node: order) {
			int max = -1;
			for (int child: graph.get(node)) if (child != info.parent[node]) {
				if (max == -1 || size[child] > size[max]) {
					max = child;
				}
			}
			for (int child: graph.get(node)) if (child != info.parent[node]) {
				ArrayList<Integer> chain;
				if (child != max) {
					chain = new ArrayList<>();
					info.chain[child] = info.chains.size();
					info.chains.add(chain);
				} else {
					int index = info.chain[node];
					chain = info.chains.get(index);
					info.chain[child] = index;
				}
				info.index[child] = chain.size();
				chain.add(child);
			}
		}

		return info;
	}
}
