import edu.princeton.cs.algs4.Picture;

public class SeamCarver {

    /**
     * create a seam carver object based on the given picture
     * 
     * @param picture
     */
    public SeamCarver(Picture picture) {

    }

    /**
     * current picture
     * 
     * @return
     */
    public Picture picture() {
        throw new RuntimeException();
    }

    /**
     * width of current picture
     * 
     * @return
     */
    public int width() {
        throw new RuntimeException();
    }

    /**
     * height of current picture
     * 
     * @return
     */
    public int height() {
        throw new RuntimeException();
    }

    /**
     * energy of pixel at column x and row y
     * 
     * @param x
     * @param y
     * @return
     */
    public double energy(int x, int y) {
        throw new RuntimeException();
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
        throw new RuntimeException();
    }

    /**
     * remove vertical seam from current picture
     * 
     * @param seam
     */
    public void removeVerticalSeam(int[] seam) {
        throw new RuntimeException();
    }
}
