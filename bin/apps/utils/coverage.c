#include "SamFile.h"
#include "SamValidation.h"

int ReadIndexedBam(const char* inputFilename,
                   const char* outputFilename,
                   const char* indexFilename)
{
   // Open the input file for reading.
   SamFile samIn;
   samIn.OpenForRead(inputFilename);

   // Open the bam index file for reading.
   samIn.ReadBamIndex(indexFilename);

   // Open the output file for writing.
   SamFile samOut;
   samOut.OpenForWrite(outputFilename);

   // Read the sam header.
   SamFileHeader samHeader;
   samIn.ReadHeader(samHeader);

   // Write the sam header.
   samOut.WriteHeader(samHeader);

   SamRecord samRecord;
   SamValidationErrors samInvalidErrors;

   // Loop through each Reference.
   for(int i = -1; i < 23; i++)
   {
      int numSectionRecords = 0;
      samIn.SetReadSection(i);

      // Keep reading records until they aren't more.
      while(samIn.ReadRecord(samHeader, samRecord))
      {
         numSectionRecords++;

         // Successfully read a record from the file, so write it.
         samOut.WriteRecord(samHeader, samRecord);
      }

      std::cout << "Reference ID " << i << " has " << numSectionRecords 
                << " records" << std::endl;
   }

   std::cout << "Number of records = " << samIn.GetCurrentRecordCount() << std::endl;

   return(0);
}
