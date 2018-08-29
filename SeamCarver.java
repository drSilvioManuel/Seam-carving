
import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Queue;

public class SeamCarver {

    private double[][] energy;
    private double[][] distTo;
    private Indegree[][] edgeTo;
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

        return new int[]{};
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
        distTo = new double[width()][height()];
        edgeTo = new Indegree[width()][height()];

        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                if (x == 0) distTo[x][y] = 0;
                else distTo[x][y] = Double.POSITIVE_INFINITY;
            }
        }
        Indegree[][] indegree = new Indegree[width()][height()];
        Iterable<Indegree> topologicalY = calculateTopologicalY(indegree);

        for (Indegree v : topologicalY) {
            for (Indegree w : adjY(v, indegree))
                relaxY(v, w);
        }

        return new int[]{};
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
    private void relaxY(Indegree v, Indegree w) {
        int y2 = w.y;
        int x2 = w.x;
        double e = energy[w.x][w.y];
        if (distTo[w.x][w.y] > distTo[v.x][v.y] + e) {
            distTo[w.x][w.y] = distTo[v.x][v.y] + e;
            edgeTo[w.x][w.y] = v;
        }
    }

    private Iterable<Indegree> calculateTopologicalY(Indegree[][] indegree) {
        // indegrees of remaining vertices
        Queue<Indegree> queue = new Queue<Indegree>();
        for (int x = 0; x < width(); x++) {
            for (int y = 0; y < height(); y++) {
                int indegreeY = indegreeY(x, y);
                indegree[x][y] = new Indegree(x, y, indegreeY);
                if (indegreeY == 0) queue.enqueue(indegree[x][y]);
            }
        }

        // vertices in topological order
        Queue<Indegree> order = new Queue<Indegree>();

        while (!queue.isEmpty()) {
            Indegree v = queue.dequeue();
            order.enqueue(v);
            for (Indegree w : adjY(v, indegree)) {
                w.indegree--;
                if (w.indegree == 0) queue.enqueue(w);
            }
        }
        return order;
    }

    private int indegreeY(int x, int y) {
        if (x == 0) return 0;
        if (y == 0 || y == width()-1) return 2;
        return 3;
    }

    private Indegree[] adjY(Indegree indegree, Indegree[][] indegrees) {
        int x = indegree.x;
        int y = indegree.y;

        if (y == height()-1) return new Indegree[]{indegrees[x+1][y-1], indegrees[x+1][y]};
        else if (y == 0) return new Indegree[]{indegrees[x+1][y], indegrees[x+1][y+1]};
        return new Indegree[]{indegrees[x+1][y-1], indegrees[x+1][y], indegrees[x+1][y+1]};
    }

    private static class Indegree {
        final int x;
        final int y;
        int indegree;

        Indegree(int x0, int y0, int indegree0) {
            x = x0;
            y = y0;
            indegree = indegree0;
        }

    }
}
