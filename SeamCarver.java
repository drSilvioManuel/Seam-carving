import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;

public class SeamCarver {

    private double[][] energy;
    private double[] distTo;
    private int[] edgeTo;
    private int width;
    private int height;
    private int [][] colorMatrix;

    /**
     * create a seam carver object based on the given picture
     * 
     * @param picture
     */
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException("The picture arg does not have to be null");
        width = picture.width();
        height = picture.height();
        energy = new double[width][height];
        colorMatrix = new int[width][height];

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                colorMatrix[x][y] = picture.getRGB(x, y);
                if (x+1 < width) colorMatrix[x+1][y] = picture.getRGB(x+1, y);
                if (y+1 < height) colorMatrix[x][y+1] = picture.getRGB(x, y+1);
                energy(x, y);
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
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++) {
                pic.setRGB(x, y, colorMatrix[x][y]);
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
     * energy of pixel at column x and row y
     * 
     * @param x
     * @param y
     * @return
     */
    public double energy(int x, int y) {
        if (x < 0 || y < 0 || x > width()-1 || y > height()-1) throw new IllegalArgumentException("Wrong x y coordinates");
        int x1 = x - 1;
        int x2 = x + 1;
        int y1 = y - 1;
        int y2 = y + 1;

        if (x == 0 || y == 0 || x == width()-1 || y == height()-1) energy[x][y] = 1000;
        else {
            int colorX1 = colorMatrix[x1][y];
            int colorX2 = colorMatrix[x2][y];
            int colorY1 = colorMatrix[x][y1];
            int colorY2 = colorMatrix[x][y2];

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

            energy[x][y] = Math.sqrt(dx + dy);
        }
        return energy[x][y];
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

        findVerticalSeam();

        energy = transposeMatrix(energy);
        colorMatrix = transposeMatrix(colorMatrix);

        flipWH();

        return edgeTo;
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
        int size = width()*height();
        distTo = new double[size+2];
        edgeTo = new int[height()];

        for (int v = 0; v < size; v++)
            distTo[v] = Double.POSITIVE_INFINITY;
        distTo[size] = 0;

        Iterable<Integer> topologicalY = calculateTopologicalY(size);

        for (int v : topologicalY) {
            for (int w : adjY(v))
                relaxY(v, w);
        }

        return edgeTo;
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
            System.arraycopy(energy[x], y+1, energy[x], y, width - y - 1);
            System.arraycopy(colorMatrix[x], y+1, colorMatrix[x], y, width - y - 1);
            energy[x][width-1] = 0;
            x++;
        }
        x = 0;
        for (int y : seam) {
            energy(x, y);
            x++;
        }
        width--;
    }

    // relax edge e
    private void relaxY(int v, int w) {
        int y2 = w / width();
        int x2 = w % width();
        double e;
        if (w >= height()*width()) {
            e = 1000;
        } else e = energy[x2][y2];
        if (distTo[w] > distTo[v] + e) {
            distTo[w] = distTo[v] + e;
            edgeTo[y2] = x2;
        }
    }

    private Iterable<Integer> calculateTopologicalY(int size) {
        // indegrees of remaining vertices
        int[] indegree = new int[size+2];
        for (int v = 0; v < size; v++) {
            indegree[v] = indegreeY(v);
        }


        // vertices in topological order
        Queue<Integer> order = new Queue<Integer>();

        // initialize queue to contain all vertices with indegree = 0
        Queue<Integer> queue = new Queue<Integer>();
        queue.enqueue(-1);

        while (!queue.isEmpty()) {
            int v = queue.dequeue();
            order.enqueue(v);
            for (int w : adjY(v)) {
                indegree[w]--;
                if (indegree[w] == 0) queue.enqueue(w);
            }
        }
        return order;
    }

    private int indegreeY(int node) {
        int y = node / width();
        int x = node % width();
        if (node == width()*height()) return 0;
        if (node == width()*height()+1) return width();
        if (y == 0) return 1;
        if (x == 0 || x == width()-1) return 2;
        return 3;
    }

    private int[] adjY(int node) {
        int y = node / width();
        int x = node % width();
        if (node == width()*height()) {
            int[] res = new int[width()];
            for (int i = 0; i < width(); i++) res[i] = i;
        } else if (node == width()*height()+1) return new int[]{};
        else if (y == height()-1) return new int[]{height()*width()+1};
        else if (x == 0) return new int[]{node+width(), node+width()+1};
        else if (x == width()-1) return new int[]{node+width()-1, node+width()};
        return new int[]{node+width()-1, node+width(), node+width()+1};
    }
}
