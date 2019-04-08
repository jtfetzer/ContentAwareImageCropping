package org.coursera;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import edu.princeton.cs.algs4.Picture;

public class SeamCarver {

	private static final double BORDER_ENERGY = 1000;
	private final List<Double> energyMatrix;
	private List<Color> mutablePicture;
	private Picture picture;
	private int width;
	private int height;
	private boolean picUpToDate = true;

	public SeamCarver(Picture picture) { // create a seam carver object based on the given picture
		if (picture == null) {
			throw new IllegalArgumentException("Constructor cannot be passed a null value.");
		}
		energyMatrix = new ArrayList<>(picture.width() * picture.height());
		width = picture.width();
		height = picture.height();
		initializeMutablePicture(picture);
		initializeEnergyMatrix();
	}

	public <T> void print(List<T> list) { // Used for debugging
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				System.out.printf(" %8.2f", list.get(i * width + j));
			}
			System.out.println();
		}
		System.out.println();
	}

	private void initializeMutablePicture(Picture picture) { // Creates a copy of the picture, so as not to mutate the input.
		mutablePicture = new ArrayList<>(picture.width() * picture.height());

		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				mutablePicture.add(picture.get(col, row));
			}
		}
	}

	private void initializeEnergyMatrix() { // Determines and sets the energy of each pixel
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				if (row == 0 || row == height - 1 || col == 0 || col == width - 1) {
					energyMatrix.add(BORDER_ENERGY);
				} else {
					energyMatrix.add(energy(col, row));
				}
			}
		}
	}

	public Picture picture() { // return current state of mutablePicture

		if (picUpToDate) {
			return picture;
		}
		final Picture pic = new Picture(width, height);

		for (int row = 0; row < height; row++) { // row
			for (int col = 0; col < width; col++) {
				pic.set(col, row, mutablePicture.get(row * width + col));
			}
		}
		picUpToDate = true;
		picture = pic;
		return picture;
	}

	public int height() { // width of current picture
		return height;
	}

	public int width() { // height of current picture
		return width;
	}

	public double energy(int col, int row) { // energy of pixel at column x and row y
		if (row < 0 || row > height - 1 || col < 0 || col > width - 1) {
			throw new IllegalArgumentException("Invalid argument");
		}
		if (row == 0 || row == (height - 1)  || col == 0 || col == (width - 1)) {
			return BORDER_ENERGY;
		}

		return dualGradient(col, row);
	}

	private double dualGradient(int col, int row) {
		return Math.sqrt(xGradientSquared(col, row) + yGradientSquared(col, row));
	}

	private double yGradientSquared(int col, int row) {
		final Color y1 = mutablePicture.get(width * (row - 1) + col);
		final Color y2 = mutablePicture.get(width * (row + 1) + col);
		return Math.pow(redEnergy(y1, y2), 2) + Math.pow(greenEnergy(y1, y2), 2) + Math.pow(blueEnergy(y1, y2), 2);
	}

	private double xGradientSquared(int col, int row) {
		final Color x1 = mutablePicture.get(width * row + (col - 1));
		final Color x2 = mutablePicture.get(width * row + (col + 1));
		return Math.pow(redEnergy(x1, x2), 2) + Math.pow(greenEnergy(x1, x2), 2) + Math.pow(blueEnergy(x1, x2), 2);
	}

	private int blueEnergy(Color c1, Color c2) {
		return c1.getBlue() - c2.getBlue();
	}

	private int greenEnergy(Color c1, Color c2) {
		return c1.getGreen() - c2.getGreen();
	}

	private int redEnergy(Color c1, Color c2) {
		return c1.getRed() - c2.getRed();
	}


	public int[] findHorizontalSeam() { // sequence of indices for horizontal seam

		final double[][] shortestPathWeight = new double[height][width];
		final int[] edgeTo = new int[energyMatrix.size()];

		final int rows = shortestPathWeight.length;
		final int cols = shortestPathWeight[0].length;
		for (int row = 0; row < rows; row++) { // Initialize shortestPathWeight matrix values to Infinity
			for (int col = 0; col < cols; col++) {
				if (col == 0) {
					shortestPathWeight[row][col] = BORDER_ENERGY;
				} else {
					shortestPathWeight[row][col] = Double.POSITIVE_INFINITY;
				} // else
			} // end for
		} // end for

		for (int col = 0; col < cols - 1; col++) { // check up to size of arrays - 1
			for (int row = 0; row < rows; row++) { // check up to number of arrays
				final double pathWeight = shortestPathWeight[row][col];
				double frontierNodeWeight, edgeValue;

				if (row != 0) {
					frontierNodeWeight = shortestPathWeight[row - 1][col + 1]; // upper right vertex shortest path value
					edgeValue = energyMatrix.get(width * (row - 1) + (col + 1)); // simulated "directed-edge" value
					if (frontierNodeWeight > pathWeight + edgeValue) {
						shortestPathWeight[row - 1][col + 1] = pathWeight + edgeValue;
						edgeTo[(row - 1) * (shortestPathWeight[0].length) + (col + 1)] = row * shortestPathWeight[0].length + col;
					} // end if
				} // end if

				frontierNodeWeight = shortestPathWeight[row][col + 1]; // center right vertex shortest path value
				edgeValue = energyMatrix.get(width * row + (col + 1));
				if (frontierNodeWeight > pathWeight + edgeValue) {
					shortestPathWeight[row][col + 1] = pathWeight + edgeValue;
					edgeTo[row * cols + (col + 1)] = row * cols + col;
				} // end if

				if (row != shortestPathWeight.length - 1) {
					frontierNodeWeight = shortestPathWeight[row + 1][col + 1]; // lower right vertex shortest path value
					edgeValue = energyMatrix.get(width * (row + 1) + (col + 1));
					if (frontierNodeWeight > pathWeight + edgeValue) {
						shortestPathWeight[row + 1][col + 1] = pathWeight + edgeValue;
						edgeTo[(row + 1) * cols + (col + 1)] = row * cols + col;
					} // end if
				}
			} // end for
		} // end for

		double min = Double.POSITIVE_INFINITY;
		int minRow = 0;
		for (int row = 0; row < shortestPathWeight.length; row++) {
			if (min > shortestPathWeight[row][cols - 1]) {
				min = shortestPathWeight[row][cols - 1];
				minRow = row;
			}
		}
		return shortestHorizontalPath(edgeTo, cols, (minRow + 1) * cols - 1);
	}

	private int[] shortestHorizontalPath(int[] nodes, int pathSize, int endNode) {

		final int[] path = new int[pathSize];
		path[pathSize - 1] = endNode;

		for (int i = pathSize - 2; i >= 0; i--) {
			path[i] = nodes[path[i + 1]];
		}

		final int[] shortestPath = new int[path.length];
		for (int i = 0; i < path.length; i++) {
			shortestPath[i] = path[i] / path.length;
		}
		return shortestPath;
	}

	public int[] findVerticalSeam() { // sequence of indices for vertical seam

		final double[][] shortestPathWeight = new double[height][width];
		final int[] edgeTo = new int[energyMatrix.size()];

		for (int row = 0; row < shortestPathWeight.length; row++) { // Initialize shortestPathWeight matrix values to Infinity
			for (int col = 0; col < shortestPathWeight[0].length; col++) {
				if (row == 0) {
					shortestPathWeight[row][col] = BORDER_ENERGY;
				} else {
					shortestPathWeight[row][col] = Double.POSITIVE_INFINITY;
				} // else
			} // end for
		} // end for

		for (int row = 0; row < shortestPathWeight.length - 1; row++) { // check up to number of arrays
			for (int col = 0; col < shortestPathWeight[0].length; col++) { // check up to size of arrays - 1
				final double pathWeight = shortestPathWeight[row][col];
				double frontierNodeWeight, edgeValue;

				if (col != 0) {
					frontierNodeWeight = shortestPathWeight[row + 1][col - 1]; // lower left vertex shortest path value
					edgeValue = energyMatrix.get(width * (row + 1) + (col - 1)); // simulated "directed-edge" value
					if (frontierNodeWeight > pathWeight + edgeValue) {
						shortestPathWeight[row + 1][col - 1] = pathWeight + edgeValue;
						edgeTo[(row + 1) * (shortestPathWeight[0].length) + (col - 1)] = row * shortestPathWeight[0].length + col;
					} // end if
				} // end if

				frontierNodeWeight = shortestPathWeight[row + 1][col]; // lower center vertex shortest path value
				edgeValue = energyMatrix.get(width * (row + 1) + (col));
				if (frontierNodeWeight > pathWeight + edgeValue) {
					shortestPathWeight[row + 1][col] = pathWeight + edgeValue;
					edgeTo[(row + 1) * shortestPathWeight[0].length + col] = row * shortestPathWeight[0].length + col;
				} // end if

				if (col != shortestPathWeight[0].length - 1) {
					frontierNodeWeight = shortestPathWeight[row + 1][col + 1]; // lower right vertex shortest path value
					edgeValue = energyMatrix.get(width * (row + 1) + (col + 1));
					if (frontierNodeWeight > pathWeight + edgeValue) {
						shortestPathWeight[row + 1][col + 1] = pathWeight + edgeValue;
						edgeTo[(row + 1) * shortestPathWeight[0].length + (col + 1)] = row * shortestPathWeight[0].length + col;
					} // end if
				}
			} // end for
		} // end for

		double min = Double.POSITIVE_INFINITY;
		int minCol = 0;
		for (int col = 0; col < shortestPathWeight[0].length; col++) {
			if (min > shortestPathWeight[shortestPathWeight.length - 1][col]) {
				min = shortestPathWeight[shortestPathWeight.length - 1][col];
				minCol = col;
			}
		}
		return shortestVerticalPath(edgeTo, shortestPathWeight.length, shortestPathWeight[0].length,
				minCol + shortestPathWeight[0].length * (shortestPathWeight.length - 1));
	}

	private int[] shortestVerticalPath(int[] nodes, int pathSize, int columns, int endNode) {

		final int[] path = new int[pathSize];
		path[pathSize - 1] = endNode;

		for (int i = pathSize - 2; i >= 0; i--) {
			path[i] = nodes[path[i + 1]];
		}

		final int[] shortestPath = new int[path.length];

		for (int i = 0; i < path.length; i++) {
			shortestPath[i] = path[i] % columns;
		}
		return shortestPath;
	}

	private boolean validVertical(int[] seam) {
		if (seam == null) {
			return false;
		}
		if (width <= 1) {
			return false;
		}
		if (seam.length != height) {
			return false;
		}
		for (int i = 0; i < seam.length; i++) {
			if (seam[i] < 0 || seam[i] >= width) {
				return false;
			}
		}
		for (int i = 0; i < seam.length - 1; i++) {
			if (Math.abs(seam[i] - seam[i + 1]) > 1) {
				return false;
			}
		}
		return true;
	}

	private boolean validHorizontal(int[] seam) {
		if (seam == null) {
			return false;
		}
		if (height <= 1) {
			return false;
		}
		if (seam.length != width) {
			return false;
		}
		for (int i = 0; i < seam.length; i++) {
			if (seam[i] < 0 || seam[i] >= height) {
				return false;
			}
		}
		for (int i = 0; i < seam.length - 1; i++) {
			if (Math.abs(seam[i] - seam[i + 1]) > 1) {
				return false;
			}
		}
		return true;
	}

	public void removeVerticalSeam(int[] seam) { // Removes vertical seam from picture.
		// contains the column pixels to remove for each row.

		if (!validVertical(seam)) {
			throw new IllegalArgumentException("Invalid seam.");
		}

		int index;
		for (int i = 0; i < seam.length; i++) {
			index = i * width + seam[i];
			energyMatrix.remove(index - i);
			mutablePicture.remove(index - i);
		}
		width--;
		picUpToDate = false;
	} // end method verticalHelper

	public void removeHorizontalSeam(int[] seam) { // Rows

		if (!validHorizontal(seam)) {
			throw new IllegalArgumentException("Invalid seam.");
		}
		int index;
		final double sum = IntStream.of(seam).sum();

		if (sum <= ((height - 1) * width) / 2) {
			for (int col = 0; col < seam.length; col++) {
				index = width * seam[col] + col;
				while (index - width >= 0) {
					energyMatrix.set(index, energyMatrix.get(index - width));
					mutablePicture.set(index, mutablePicture.get(index - width));
					index -= width;
				}
			}
			for (int i = 0; i < seam.length; i++) {
				energyMatrix.remove(0);
				mutablePicture.remove(0);
			}
		} else {
			for (int col = 0; col < seam.length; col++) {
				index = width * seam[col] + col;
				while (index + width < width * height) {
					energyMatrix.set(index, energyMatrix.get(index + width));
					mutablePicture.set(index, mutablePicture.get(index + width));
					index += width;
				}
			}
			for (int i = 0; i < seam.length; i++) {
				index = height * width - 1 - i;
				energyMatrix.remove(index);
				mutablePicture.remove(index);
			}
		}
		height--;
		picUpToDate = false;
	}
}
