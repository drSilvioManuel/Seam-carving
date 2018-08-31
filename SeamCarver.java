import edu.princeton.cs.algs4.MinPQ;
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;

import java.util.HashMap;

public class SeamCarver {

    private static final int BLUE = 0xff;
    private static final int GREEN = 0xff00;
    private static final int RED = 0xff0000;
    private static final int MIN_ACCEPTABLE_CELLS = 3;
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

            int blueX1 = colorX1 & BLUE;
            int greenX1 = (colorX1 & GREEN) >> 8;
            int redX1 = (colorX1 & RED) >> 16;

            int blueX2 = colorX2 & BLUE;
            int greenX2 = (colorX2 & GREEN) >> 8;
            int redX2 = (colorX2 & RED) >> 16;

            int blueY1 = colorY1 & BLUE;
            int greenY1 = (colorY1 & GREEN) >> 8;
            int redY1 = (colorY1 & RED) >> 16;

            int blueY2 = colorY2 & BLUE;
            int greenY2 = (colorY2 & GREEN) >> 8;
            int redY2 = (colorY2 & RED) >> 16;

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

        if (cols < MIN_ACCEPTABLE_CELLS || rows < MIN_ACCEPTABLE_CELLS) return new int[]{};

        distTo = new double[rows][cols];

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                if (row == 0) distTo[row][col] = 0;
                else distTo[row][col] = Double.POSITIVE_INFINITY;
            }
        }
        Node[][] nodes = new Node[rows][cols];
        Iterable<Node> topologicalY = calculateTopologicalY(nodes);
        SemiGraph graph = new SemiGraph(rows, cols);
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

        int[][] copyColor = new int[height()][width() - 1];
        double[][] copyEnergy = new double[height()][width() - 1];

        for (int row = 0; row < height(); row++) {
            System.arraycopy(colorMatrix[row], 0, copyColor[row], 0, seam[row]);
            System.arraycopy(colorMatrix[row], seam[row] + 1, copyColor[row], seam[row], height - seam[row] - 1);

            System.arraycopy(energy[row], 0, copyEnergy[row], 0, seam[row]);
            System.arraycopy(energy[row], seam[row] + 1, copyEnergy[row], seam[row], height - seam[row] - 1);
        }
        colorMatrix = copyColor;
        energy = copyEnergy;

        for (int row = 0; row < height(); row++) energy[row][seam[row]] = energy(seam[row], row);

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
        private static final int INIT_HASH = 17;
        private static final int FOLLOW_HASH = 31;
        private final int row;
        private final int col;
        private int indegree;
        private double power;

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
            int result = INIT_HASH;
            result = FOLLOW_HASH * result + row;
            result = FOLLOW_HASH * result + col;
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

        private final HashMap<Node, MinPQ<Node>> map = new HashMap<>();
        private Node lastNode;
        private final int row;
        private final int col;

        SemiGraph(int r, int c) {
            row = r;
            col = c;
        }


        void addEdge(Node from, Node to, double weight) {
            MinPQ<Node> pq = map.get(to);
            if (pq == null) {
                pq = new MinPQ<>(3, (o1, o2) -> {
                    double res = o1.power - o2.power;
                    if (res > 0) return 1;
                    else if (res < 0) return -1;
                    else return 0;
                });
                map.put(to, pq);
            }
            from.power = weight;
            pq.insert(from);
            if (to.row == row-1) {
                if (lastNode == null) lastNode = from;
                else if (from.power < lastNode.power) lastNode = from;
            }
        }

        int[] getChain(int size) {
            int[] res = new int[size];
            int i = size;
            res[--i] = 0;
            for (Node next = lastNode; next != null && map.get(next) != null; next = map.get(next).delMin()) {
                res[--i] = next.col;
            }
            res[size-1] = res[size-2]; // last elements have same power, so assigh as prev
            res[0] = res[1]; // first elements have same power, so assigh as next
            return res;
        }
    }
}
