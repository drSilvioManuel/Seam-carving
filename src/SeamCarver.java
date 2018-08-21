import edu.princeton.cs.algs4.Picture;

public class SeamCarver {

    private final Picture picture;
    private final double[][] energy;

    /**
     * create a seam carver object based on the given picture
     * 
     * @param picture
     */
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException("The picture arg does not have to be null");
        this.picture = picture;
        this.energy = new double[picture.width()][picture.height()];
    }

    /**
     * current picture
     * 
     * @return
     */
    public Picture picture() {
        return picture;
    }

    /**
     * width of current picture
     * 
     * @return
     */
    public int width() {
        return picture.width();
    }

    /**
     * height of current picture
     * 
     * @return
     */
    public int height() {
        return picture.height();
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
            int colorX1 = picture.getRGB(x1, y);
            int colorX2 = picture.getRGB(x2, y);
            int colorY1 = picture.getRGB(x, y1);
            int colorY2 = picture.getRGB(x, y2);

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
        throw new RuntimeException();
    }

    /**
     * sequence of indices for vertical seam
     * 
     * @return
     */
    public int[] findVerticalSeam() {
        throw new RuntimeException();
    }

    /**
     * remove horizontal seam from current picture
     * 
     * @param seam
     */
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("The horizontal seam does not have to be null");
        throw new RuntimeException();
    }

    /**
     * remove vertical seam from current picture
     * 
     * @param seam
     */
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("The vertical seam does not have to be null");
        throw new RuntimeException();
    }
}
