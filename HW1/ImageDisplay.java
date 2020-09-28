import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;

import javax.swing.*;


public class ImageDisplay {

	public JFrame frame;
	public JLabel lbIm1;
	public BufferedImage  imgOne;
	public static int width = 512;
	public static int height = 512;
	public boolean doScale = false;
	public boolean uniQuanti = false;
	public boolean logQuanti = false;
	public ArrayList<Integer> pixelsList = new ArrayList<>();
	public ArrayList<Integer> UnipixelsList = new ArrayList<>();
	public ArrayList<Integer> rList = new ArrayList<>();
	public ArrayList<Integer> gList = new ArrayList<>();
	public ArrayList<Integer> bList = new ArrayList<>();
	public int[] uniBuckets = null;
	public int[] logBuckets = null;

	public void showIms(String[] args){

		// Read in the specified image
		
		try
		{
			int frameLength = width*height*3;
			String imgPath = args[0];
			File file = new File(imgPath);
			RandomAccessFile raf = new RandomAccessFile(file, "r");
			raf.seek(0);
	
			long len = frameLength;
			byte[] bytes = new byte[(int) len];

			raf.read(bytes);
			
//			Calculate the original Pixels
	        int ind = 0;
			for(int y = 0; y < height; y++)
			{
				for(int x = 0; x < width; x++)
				{
					byte a = 0;
					byte r = bytes[ind];
					byte g = bytes[ind+height*width];
					byte b = bytes[ind+height*width*2]; 
					
					int R = r & 0xff;
					rList.add(R);
					int G = g & 0xff;
					gList.add(G);
					int B = b & 0xff;
					bList.add(B);

					int pix = 0xff000000 | ((R) << 16) | ((G) << 8) | (B);
					pixelsList.add(pix);
					
					ind++;

				}
			}
			
//			transfer pixelsList to array pixels[]
			int[] pixels=new int[pixelsList.size()];  
	        for(int i=0;i<pixelsList.size();i++){  
	        	pixels[i] = pixelsList.get(i);  
	        }  
	        
//			Filter --- Median Filtering
			
			int[] newPixels = medianFiltering(pixels, width, height);
//			update original pixel list and RGB list
			pixelsList.clear();
			for(int i = 0; i < newPixels.length; i++) {
				pixelsList.add(newPixels[i]);
			}
			calculateRGB(width, height, newPixels);
			
//			Check Do Scale Or Not
	        Double scalNum = Double.valueOf(args[1]); 
	        if(scalNum > 1.0 || scalNum <= 0) {
	        	System.out.println("Invalid scale number. Scale Value should be greater than 0.0 but less than or equal to 1.0. ");
	        }
	        if(scalNum == 1.0) {
	        	doScale = false;
			}else {
				doScale = true;

			}

//			Check Do Uniform Quantization Or Not
	        int QuantiPivot = Integer.valueOf(args[3]);
	        int uniQuantiNum = Integer.valueOf(args[2]);
//	        
	        if(uniQuantiNum > 8 || uniQuantiNum < 1) {
	        	System.out.println("Invalid Quantization number. Quantization Value should be greater or equal to 1 and less than or equal to 8. ");
	        }
	        
	        if(QuantiPivot < -1 || QuantiPivot > 255) {
	        	System.out.println("Invalid Quantization mode. Mode Value should be greater or equal to -1 and less than or equal to 255. ");
	        }
	        if(uniQuantiNum == 8) {
	        	uniQuanti = false;
	        	logQuanti = false;
	        }else {
	        	double valPerChannel = Math.pow(2, uniQuantiNum);
 				double slotSize = 256 / valPerChannel;
 				if(QuantiPivot == -1) {
// 					Do Uniform Quantization
 					uniQuanti = true;
 					uniBuckets = getBucketsForUniformQuantization(slotSize);
 				}else {
// 					Do Log Quantization
 					logQuanti = true;
 					uniBuckets = getBucketsForUniformQuantization(slotSize);
 					logBuckets = getBucketsForLogQuantization(QuantiPivot, uniBuckets);
 				}
	        }
//

	        if(doScale == true) {
	        	
	        	int newWidth = new Double(width * scalNum).intValue();
				int newHeight = new Double(height * scalNum).intValue();
				imgOne = new BufferedImage(newWidth, newHeight, BufferedImage.TYPE_INT_RGB);
		        int[] newnewPixels = resizeBilinear(newPixels, width, height, newWidth, newHeight);
		        calculateRGB(newWidth, newHeight, newnewPixels);
//		        display(newPixels, imgOne, newWidth, newHeight);
		        if(logQuanti == true ) {
		        	doQuantization(newWidth, newHeight, logBuckets);
//	        		transfer UnipixelsList to array Unipixels[]
	        		int[] Unipixels=new int[UnipixelsList.size()];  
	                for(int i=0; i < UnipixelsList.size(); i++){ 
	                	Unipixels[i] = UnipixelsList.get(i);  
	                }  
		    		display(Unipixels, imgOne, newHeight, newHeight);
		        }
		        else if(uniQuanti == true) {
	        		doQuantization(newWidth, newHeight, uniBuckets);
//	        		transfer UnipixelsList to array Unipixels[]
	        		int[] Unipixels=new int[UnipixelsList.size()];  
	                for(int i=0; i < UnipixelsList.size(); i++){ 
	                	Unipixels[i] = UnipixelsList.get(i);  
	                }  
		    		display(Unipixels, imgOne, newHeight, newHeight);
		        }else {
		        	display(newnewPixels, imgOne, newHeight, newHeight);
		        }
	        }else {
	        	imgOne = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
				if(logQuanti == true) {
					doQuantization(width, height, logBuckets);
//	        		transfer UnipixelsList to array Unipixels[]
	        		int[] Unipixels=new int[UnipixelsList.size()];  
	                for(int i=0; i < UnipixelsList.size(); i++){ 
	                	Unipixels[i] = UnipixelsList.get(i);  
	                }  
		    		display(Unipixels, imgOne, width, height);	        	
				}
		        else if(uniQuanti == true) {
	        		doQuantization(width, height, uniBuckets);
//	        		transfer UnipixelsList to array Unipixels[]
	        		int[] Unipixels=new int[UnipixelsList.size()];  
	                for(int i=0; i < UnipixelsList.size(); i++){ 
	                	Unipixels[i] = UnipixelsList.get(i);  
	                }  
		    		display(Unipixels, imgOne, width, height);
		        }else {
		        	display(newPixels, imgOne, width, height);
		        }
	        }

		}catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
		
		
		
	}
	
	/**
	 * @param newWidth		image width
	 * @param newHeight		image height
	 * @param buckets		quantization intervals 
	 */
	public void doQuantization(int newWidth, int newHeight, int[] buckets) {
		int index = 0;
		for(int y = 0; y < newHeight; y++)
		{
			for(int x = 0; x < newWidth; x++)
			{	
				int[] quantizedRGB = quantize(rList.get(index), gList.get(index), bList.get(index), buckets);
				int newR = quantizedRGB[0];
				int newG = quantizedRGB[1];
				int newB = quantizedRGB[2];
				int uniQuantiPixel = 0xff000000 | ((newR) << 16) | ((newG) << 8) | (newB);
				UnipixelsList.add(uniQuantiPixel);
//				imgOne.setRGB(x, y, uniQuantiPixel);
				index++;
			}
		}

	}
	
	/**	calculate RGB value image pixels
	 * @param currentWidth		image width
	 * @param currentHeight		image height 
	 * @param currentPixel		image pixel
	 */
	public void calculateRGB(int currentWidth, int currentHeight, int[] currentPixel) {
		int ind = 0;
		gList.clear();
		rList.clear();
		bList.clear();
		for(int y = 0; y < currentHeight; y++)
		{
			for(int x = 0; x < currentWidth; x++)
			{
				int pixel = currentPixel[ind]; 
				int r = (pixel & 0xff0000) >> 16;
				int g = (pixel & 0xff00) >> 8;
				int b = (pixel & 0xff);
//				int R = r & 0xff;
				rList.add(r);
//				int G = g & 0xff;
				gList.add(g);
//				int B = b & 0xff;
				bList.add(b);

				ind++;

			}
		}
	}
	
	
	
	/**  Show pictures 
	 * @param pixels		image pixels
	 * @param imageShow		image ready to show 
	 * @param width			image width
	 * @param height		image height
	 */
	public void display(int[] pixels, BufferedImage imageShow, int width, int height) {
		
		int ind = 0;
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				imageShow.setRGB(x,y,pixels[ind]);
				ind++;
			}
		}
		
		// Use label to display the image
		frame = new JFrame();
		GridBagLayout gLayout = new GridBagLayout();
		frame.getContentPane().setLayout(gLayout);

		lbIm1 = new JLabel(new ImageIcon(imageShow));

		GridBagConstraints c = new GridBagConstraints();
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		c.weightx = 0.5;
		c.gridx = 0;
		c.gridy = 0;

		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridx = 0;
		c.gridy = 1;
		frame.getContentPane().add(lbIm1, c);

		frame.pack();
		frame.setVisible(true);
	}


	public static void main(String[] args) {
		ImageDisplay ren = new ImageDisplay();
		ren.showIms(args);
		
	}
	
	
	/**
	 * @param 	Divide intervals  
	 * @return    uniform buckets for center bit per channel
	 */
	private int[] getBucketsForUniformQuantization(double interval) {
		LinkedList<Integer> intervalList = new LinkedList<Integer>();
		double trueValue = 0;
		int intValue = 0;
		intervalList.add(intValue);
		while(intValue <= 255){
			trueValue = trueValue + interval;
			intValue = (int) Math.round(trueValue);
			if(intValue > 255){
				break;
			}
			intervalList.add(intValue);
		}
		intervalList.add(255);
		int[] buckets = new int[intervalList.size()];
		for(int i = 0; i < intervalList.size(); i++) {
			buckets[i] = intervalList.get(i);
		}
		return buckets;
	}
	
	/**
	 * @param pivot    pivot for log quantization 
	 * @param uniBuckets	uniform buckets for center bit per channel
	 * @return 		log buckets 
	 */
	private int[] getBucketsForLogQuantization(int pivot, int[] uniBuckets) {
		
		LinkedList<Integer> logInterval = new LinkedList<Integer>();
//		get log value for uniform buckets 
//		start with i = 1 because uniBuckets[0] is always equal to 0

		LinkedList<Double> logVals = new LinkedList<Double>();
		for(int i = 1; i < uniBuckets.length; i++) {
			double logVal =  Math.log(uniBuckets[i]);
			logVals.add(logVal);
		}
		logVals.add(Math.log(256));

//		get each block of log scale 
		double log256 = logVals.get(logVals.size() - 1);
		for(int i = logVals.size() - 2; i >= 0; i--) {
			double midVal = logVals.get(i) / log256;
			int midInterval = (int) Math.round( (1 - midVal) * 256);
			logInterval.add(midInterval);

		}
		logInterval.add(255);
//		scale them to match 256 centered at pivot
		int intervalNum = 0;
		HashSet<Integer> actualInterval = new HashSet<Integer>(); 
		actualInterval.add(0);
		for(int i = 0; i < logInterval.size() - 1; i++) {
			if(intervalNum <= uniBuckets.length - 1) {
				int largeVal = pivot + logInterval.get(i);
				if(largeVal <= 255) {
//					System.out.println(largeVal);
					actualInterval.add(largeVal);
					intervalNum++;
					if(intervalNum >= uniBuckets.length - 1 ) {
						break;
					}
				}
				int smallerVal = pivot - logInterval.get(i);
				if(smallerVal >= 1) {
					actualInterval.add(smallerVal);
					intervalNum++;
					if(intervalNum >= uniBuckets.length - 1) {
						break;
					}
				}
			}
		}
		
		actualInterval.add(255);
		int[] buckets = new int[actualInterval.size()];
//		System.out.println(actualInterval);
		int index = 0;
		for(int interval : actualInterval) {
			buckets[index] = interval;
			index++;
		}
		Arrays.sort(buckets);
		return buckets;

	}
	
	/**
	 * @param R  	red value
	 * @param G		green value
	 * @param B		blue value
	 * @param buckets	 quantization interval 
	 * @return
	 */
	private int[] quantize(int rVal, int gVal, int bVal, int[] buckets) {
		
		for(int i=0; i < buckets.length-1; i++) {
			if(rVal >= buckets[i] && rVal <= buckets[i+1]){				
				int meanVal = (int) Math.round((buckets[i] + buckets[i+1])/2.0);
				if(rVal < meanVal){
					rVal = buckets[i];
				}else{
					rVal = buckets[i+1];
				}
				break;
			}
		}
//		corner case 
		if(rVal > 255){
			rVal = 255;
		}else if(rVal < 0){
			rVal = 0;
		}

		for(int i=0; i < buckets.length-1; i++){
			if(gVal >= buckets[i] && gVal <= buckets[i+1]){				
				int meanVal = (int) Math.round((buckets[i] + buckets[i+1])/2.0);
				if(gVal < meanVal){
					gVal = buckets[i];
				}else{
					gVal = buckets[i+1];
				}
				break;
			}
		}
		if(gVal > 255){
			gVal = 255;
		}else if(gVal < 0){
			gVal = 0;
		}
		
		for(int i=0; i < buckets.length-1; i++){
			if(bVal >= buckets[i] && bVal <= buckets[i+1]){				
				int meanVal = (int) Math.round((buckets[i] + buckets[i+1])/2.0);
				if(bVal < meanVal){
					bVal = buckets[i];
				}else{
					bVal = buckets[i+1];
				}
				break;
			}
		}
		if(bVal > 255){
			bVal = 255;
		}else if(bVal < 0){
			bVal = 0;
		}

		int[] quantizedRGB = new int[]{rVal, gVal, bVal};
		return quantizedRGB;
	}
	
	/**
	 * @param pixels get pixels of current pictures (after filtering)
	 * @param w  image original width
	 * @param h  image original height
	 * @param w2  image scaled width
	 * @param h2  image scaled height 
	 * @return
	 */
	public int[] resizeBilinear(int[] pixels, int w, int h, int w2, int h2) {
	    int[] temp = new int[w2*h2] ;
	    int a, b, c, d, x, y, index ;
	    float x_ratio = ((float)(w-1)) / w2 ;
	    float y_ratio = ((float)(h-1)) / h2 ;
	    float x_diff, y_diff, blue, red, green ;
	    int offset = 0 ;
	    for (int i = 0; i < h2; i++) {
	        for (int j = 0; j < w2; j++) {
	            x = (int)(x_ratio * j) ;
	            y = (int)(y_ratio * i) ;
	            x_diff = (x_ratio * j) - x ;
	            y_diff = (y_ratio * i) - y ;
	            index = (y * w+x) ;                
	            a = pixels[index] ;
	            b = pixels[index+1] ;
	            c = pixels[index+w] ;
	            d = pixels[index+w+1] ;

	            // Yb = Ab(1-w)(1-h) + Bb(w)(1-h) + Cb(h)(1-w) + Db(wh)
	            blue = (a&0xff)*(1-x_diff)*(1-y_diff) + (b&0xff)*(x_diff)*(1-y_diff) +
	                   (c&0xff)*(y_diff)*(1-x_diff)   + (d&0xff)*(x_diff*y_diff);

	            green = ((a>>8)&0xff)*(1-x_diff)*(1-y_diff) + ((b>>8)&0xff)*(x_diff)*(1-y_diff) +
	                    ((c>>8)&0xff)*(y_diff)*(1-x_diff)   + ((d>>8)&0xff)*(x_diff*y_diff);

	            red = ((a>>16)&0xff)*(1-x_diff)*(1-y_diff) + ((b>>16)&0xff)*(x_diff)*(1-y_diff) +
	                  ((c>>16)&0xff)*(y_diff)*(1-x_diff)   + ((d>>16)&0xff)*(x_diff*y_diff);

	            temp[offset++] = 
	                    0xff000000 | 
	                    ((((int)red)<<16)&0xff0000) |
	                    ((((int)green)<<8)&0xff00) |
	                    ((int)blue) ;
	        }
	    }
	    return temp ;
	}

	
	/**
	 * @param pixel  get original pixels of the original image 
	 * @param w   get image width
	 * @param h   get image height 
	 * @return   return image pixels after median filtering 
	 */
	private int[] medianFiltering(int[] pixel, int w, int h) {
        int[] newPixel = new int[w * h];
        int[] tempR = new int[9];
        int[] tempG = new int[9];
        int[] tempB = new int[9];
//        ColorModel cm = ColorModel.getRGBdefault();
        int r;
        int g;
        int b;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if (x == 0 || x == w - 1 || y == 0 || y == h - 1) {
                    newPixel[y * w + x] = pixel[y * w + x];
                    continue;
                }
//              compute red value 
                tempR[0] = rList.get(x - 1 + (y - 1) * w);
                tempR[1] = rList.get(x + (y - 1) * w);
                tempR[2] = rList.get(x + 1 + (y - 1) * w);
                tempR[3] = rList.get(x - 1 + y * w);
                tempR[4] = rList.get(x + y * w);
                tempR[5] = rList.get(x + 1 + y * w);
                tempR[6] = rList.get(x - 1 + (y + 1) * w);
                tempR[7] = rList.get(x + (y + 1) * w);
                tempR[8] = rList.get(x + 1 + (y + 1) * w);
                r = getMedianValue(tempR);

//              compute green value
                tempG[0] = gList.get(x - 1 + (y - 1) * w);
                tempG[1] = gList.get(x + (y - 1) * w);
                tempG[2] = gList.get(x + 1 + (y - 1) * w);
                tempG[3] = gList.get(x - 1 + y * w);
                tempG[4] = gList.get(x + y * w);
                tempG[5] = gList.get(x + 1 + y * w);
                tempG[6] = gList.get(x - 1 + (y + 1) * w);
                tempG[7] = gList.get(x + (y + 1) * w);
                tempG[8] = gList.get(x + 1 + (y + 1) * w);
                g = getMedianValue(tempG);
                
//              compute blue value
                tempB[0] = bList.get(x - 1 + (y - 1) * w);
                tempB[1] = bList.get(x + (y - 1) * w);
                tempB[2] = bList.get(x + 1 + (y - 1) * w);
                tempB[3] = bList.get(x - 1 + y * w);
                tempB[4] = bList.get(x + y * w);
                tempB[5] = bList.get(x + 1 + y * w);
                tempB[6] = bList.get(x - 1 + (y + 1) * w);
                tempB[7] = bList.get(x + (y + 1) * w);
                tempB[8] = bList.get(x + 1 + (y + 1) * w);
                b = getMedianValue(tempB);
                
                newPixel[y * w + x] = 255 << 24 | r << 16 | g << 8 | b;
            }
        }
        return newPixel;
    }
	
	/**  Auxiliary Function of function :  int[] medianFiltering(int[] pixel, int w, int h)
	 * @param array  array that to be sorted 
	 * @return  median of this array
	 */
	public int getMedianValue(int[] array) {
		int median = 0;
		Arrays.sort(array);
		int len = array.length;
	    if(len % 2 == 0){
	    	median = array[((len - 1)/2)];
	    }else {
	    	median = new Double((array[len/2 -1] + array[len/2])/ 2.0).intValue();
	    }
		return median;
	}
}
