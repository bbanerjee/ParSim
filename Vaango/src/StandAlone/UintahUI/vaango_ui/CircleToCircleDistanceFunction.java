package vaango_ui;

public class CircleToCircleDistanceFunction implements DistanceFunction {
  
  // The first three coordinates are x, y, z and the last is radius
  @Override
  public double distance(double[] p1, double[] p2) {
      double d = 0;

      for (int i = 0; i < 2; i++) {
          double diff = (p1[i] - p2[i]);
          d += diff * diff;
      }
      d = Math.sqrt(d) - (p1[3]+p2[3]);
      if (d < 0.0) d = 0.001*p1[3];
      return d*d;
  }
  
  //The first three coordinates are x, y, z and the last is radius
  @Override
  public double distanceToRect(double[] point, double[] min, double[] max) {
      double d = 0;

      for (int i = 0; i < 2; i++) {
          double diff = 0;
          if (point[i]+point[3] > max[i]) {
              diff = (point[i]+point[3] - max[i]);
          }
          else if (point[i]-point[3] < min[i]) {
              diff = (point[i]-point[3] - min[i]);
          }
          d += diff * diff;
      }

      return d;
  }
}