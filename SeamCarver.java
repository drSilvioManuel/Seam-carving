import edu.princeton.cs.algs4.MinPQ;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;

import java.util.HashMap;
import java.util.Map;

public class SeamCarver {

    private int [][] colorMatrix;
    private double[][] energy;
    private double[][] distTo;
    private int width;
    private int height;

    /**
     * create a seam carver object based on the given picture
     * 
     * @param picture
     */
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException("The picture arg does not have to be null");
        width = picture.width();
        height = picture.height();
        energy = new double[height][width];
        colorMatrix = new int[height][width];

        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                colorMatrix[row][col] = picture.getRGB(col, row);
                if (col+1 < width) colorMatrix[row][col+1] = picture.getRGB(col+1, row);
                if (row+1 < height) colorMatrix[row+1][col] = picture.getRGB(col, row+1);
                energy[row][col] = energy(col, row);
            }
        }
    }

    /**
     * current picture
     * 
     * @return
     */
    public Picture picture() {
        Picture pic = new Picture(width, height);
        for (int row = 0; row < height; row++)
            for (int col = 0; col < width; col++) {
                pic.setRGB(col, row, colorMatrix[row][col]);
            }
        return pic;
    }

    /**
     * width of current picture
     * 
     * @return
     */
    public int width() {
        return width;
    }

    /**
     * height of current picture
     * 
     * @return
     */
    public int height() {
        return height;
    }

    /**
     * energy of pixel at column row and row col
     * 
     * @param x
     * @param y
     * @return
     */
    public double energy(int x, int y) {
        if (x < 0 || y < 0 || x > width()-1 || y > height()-1) throw new IllegalArgumentException("Wrong row col coordinates");

        if (x == 0 || y == 0 || x == width()-1 || y == height()-1) return 1000;
        else {

            int x1 = x - 1;
            int x2 = x + 1;
            int y1 = y - 1;
            int y2 = y + 1;

            int colorX1 = colorMatrix[y][x1];
            int colorX2 = colorMatrix[y][x2];
            int colorY1 = colorMatrix[y1][x];
            int colorY2 = colorMatrix[y2][x];

            int blueX1 = colorX1 & 0xff;
            int greenX1 = (colorX1 & 0xff00) >> 8;
            int redX1 = (colorX1 & 0xff0000) >> 16;

            int blueX2 = colorX2 & 0xff;
            int greenX2 = (colorX2 & 0xff00) >> 8;
            int redX2 = (colorX2 & 0xff0000) >> 16;

            int blueY1 = colorY1 & 0xff;
            int greenY1 = (colorY1 & 0xff00) >> 8;
            int redY1 = (colorY1 & 0xff0000) >> 16;

            int blueY2 = colorY2 & 0xff;
            int greenY2 = (colorY2 & 0xff00) >> 8;
            int redY2 = (colorY2 & 0xff0000) >> 16;

            int bx = blueX2 - blueX1;
            int gx = greenX2 - greenX1;
            int rx = redX2 - redX1;

            int by = blueY2 - blueY1;
            int gy = greenY2 - greenY1;
            int ry = redY2 - redY1;

            double dx = bx*bx + gx*gx + rx*rx;
            double dy = by*by + gy*gy + ry*ry;

            return Math.sqrt(dx + dy);
        }
    }

    /**
     * sequence of indices for horizontal seam
     * 
     * @return
     */
    public int[] findHorizontalSeam() {

        energy = transposeMatrix(energy);
        colorMatrix = transposeMatrix(colorMatrix);

        flipWH();

        int[] res = findVerticalSeam();

        energy = transposeMatrix(energy);
        colorMatrix = transposeMatrix(colorMatrix);

        flipWH();

        return res;
    }

    private static double[][] transposeMatrix(double[][] m) {
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    private static int[][] transposeMatrix(int[][] m) {
        int[][] temp = new int[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    /**
     * sequence of indices for vertical seam
     *
     * @return
     */
    public int[] findVerticalSeam() {
        int rows = height();
        int cols = width();
        distTo = new double[rows][cols];

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                if (row == 0) distTo[row][col] = 0;
                else distTo[row][col] = Double.POSITIVE_INFINITY;
            }
        }
        Node[][] nodes = new Node[rows][cols];
        Iterable<Node> topologicalY = calculateTopologicalY(nodes);
        SemiGraph graph = new SemiGraph();
        for (Node v : topologicalY) {
            for (Node w : adjY(v, nodes))
                relaxY(v, w, graph);
        }

        return graph.getChain(rows);
    }

    /**
     * remove horizontal seam from current picture
     *
     * @param seam
     */
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("The horizontal seam does not have to be null");
        energy = transposeMatrix(energy);
        colorMatrix = transposeMatrix(colorMatrix);

        flipWH();

        removeVerticalSeam(seam);

        energy = transposeMatrix(energy);
        colorMatrix = transposeMatrix(colorMatrix);

        flipWH();
    }

    private void flipWH() {
        width = height + width;
        height = width - height;
        width = width - height;
    }

    /**
     * remove vertical seam from current picture
     *
     * @param seam
     */
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("The vertical seam does not have to be null");
        int x = 0;
        for (int y : seam) {
            System.arraycopy(energy[y], x+1, energy[y], x, width - x - 1);
            System.arraycopy(colorMatrix[y], x+1, colorMatrix[y], x, width - x - 1);
            energy[x][width-1] = 0;
            x++;
        }
        x = 0;
        for (int y : seam) {
            energy[y][x] = energy(x, y);
            x++;
        }
        width--;
    }

    // relax edge e
    private void relaxY(Node v, Node w, SemiGraph graph) {

        double e = energy[w.row][w.col];
        if (distTo[w.row][w.col] > distTo[v.row][v.col] + e) {
            distTo[w.row][w.col] = distTo[v.row][v.col] + e;
            graph.addEdge(v, w, distTo[v.row][v.col] + e);
        }
    }

    private Iterable<Node> calculateTopologicalY(Node[][] nodes) {
        // indegrees of remaining vertices
        Queue<Node> queue = new Queue<Node>();
        for (int row = 0; row < height(); row++) {
            for (int col = 0; col < width(); col++) {
                int indegreeY = indegreeY(row, col);
                nodes[row][col] = new Node(row, col, indegreeY);
                if (indegreeY == 0) queue.enqueue(nodes[row][col]);
            }
        }

        // vertices in topological order
        Queue<Node> order = new Queue<Node>();

        while (!queue.isEmpty()) {
            Node v = queue.dequeue();
            order.enqueue(v);
            for (Node w : adjY(v, nodes)) {
                w.indegree--;
                if (w.indegree == 0) queue.enqueue(w);
            }
        }
        return order;
    }

    private int indegreeY(int row, int col) {
        if (row == 0) return 0;
        else if (col == 0 || col == width()-1) return 2;
        else return 3;
    }

    private Node[] adjY(Node node, Node[][] nodes) {
        int row = node.row;
        int col = node.col;
        if (row == height()-1) return new Node[]{};
        else if (col == width()-1) return new Node[]{nodes[row+1][col-1], nodes[row+1][col]};
        else if (col == 0) return new Node[]{nodes[row+1][col], nodes[row+1][col+1]};
        return new Node[]{nodes[row+1][col-1], nodes[row+1][col], nodes[row+1][col+1]};
    }

    private static class Node {
        final int row;
        final int col;
        int indegree;
        Node from;

        Node(int r, int c, int indegree0) {
            row = r;
            col = c;
            indegree = indegree0;
        }

        @Override
        public String toString() {
            return String.format("Node(row %d, col %d, indegree %d)", row, col, indegree);
        }

        @Override
        public int hashCode() {
            int result = 17;
            result = 31 * result + row;
            result = 31 * result + col;
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            boolean result = false;
            if (obj == null || obj.getClass() != getClass()) {
                result = false;
            } else {
                Node person = (Node) obj;
                if (this.row == person.row && this.col == person.col) {
                    result = true;
                }
            }
            return result;
        }
    }

    private static class SemiGraph {

        Map<Node, MinPQ<Node>> map = new HashMap<>();
        MinPQ<Node> lastNodes;


        void addEdge(Node from, Node to, double weight) {
            MinPQ<Node> pq = map.get(to);
            if (pq == null) {
                pq = new MinPQ<>((o1, o2) -> o1.indegree - o2.indegree);
                map.put(to, pq);
            }
            from.indegree = (int) weight;
            from.from = to;
            pq.insert(from);
            lastNodes = pq;
        }

        int[] getChain(int size) {
            int[] res = new int[size];
            Node root = lastNodes.delMin();
            int i = 0;
            res[i++] = root.from.col;
            res[i++] = root.col;
            Node next;
            while ((next = map.get(root).delMin()) != null) {
                res[i++] = next.col;
                root = next;
                if (map.get(root) == null) break;
            }
            return res;
        }
    }
}
