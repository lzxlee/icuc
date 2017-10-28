import com.mathworks.toolbox.javabuilder.*;
import makesqr.*;

class wavreq
{
   public static void main(String[] args)
   {
      MWNumericArray n = null;
      Object[] result = null;
      Class1 theMagic = null;

      if (args.length == 0)
      {
        System.out.println("Error: must input a positive integer");
        return;
      }

      try
      {
         n = new MWNumericArray(Double.valueOf(args[0]),
                                      MWClassID.DOUBLE);

         theMagic = new Class1();

         result = theMagic.makesqr(1, n);
         System.out.println(result[0]);
      }
      catch (Exception e)
      {
         System.out.println("Exception: " + e.toString());
      }
      finally
      {
         MWArray.disposeArray(n);
         MWArray.disposeArray(result);
         theMagic.dispose();
      }
   }
}