//#ifndef __SAM_FILE_H__
//#define __SAM_FILE_H__
//
//#include "SamStatus.h"
//#include "InputFile.h"
//#include "SamFileHeader.h"
//#include "SamRecord.h"
//#include "GenericSamInterface.h"
//#include "BamIndex.h"
//#include "SamStatistics.h"
//
//class SamFile
//{
//public:
//    enum OpenType {READ, WRITE};
//
//    /// Enum for indicating the type of sort for the file.
//    enum SortedType {
//        UNSORTED = 0, ///< file is not sorted.
//        FLAG,         ///< SO flag from the header indicates the sort type.
//        COORDINATE,   ///< file is sorted by coordinate.
//        QUERY_NAME    ///< file is sorted by queryname.
//    };
//    
//    /// Default Constructor.
//    SamFile();
//
//    /// Constructor that opens the specified file based on the specified mode
//    /// (READ/WRITE).  Default is READ.
//    /// \param filename name of the file to open.
//    /// \param mode mode to use for opening the file (defaults to READ).
//    SamFile(const char* filename, OpenType mode = READ);
//
//    virtual ~SamFile();
//   
//    /// Open a sam/bam file for reading with the specified filename.
//    /// \param  filename: the sam/bam file to open for reading.
//    /// \return true = success; false = failure.   
//    bool OpenForRead(const char * filename);
//
//    /// Open a sam/bam file for writing with the specified filename.
//    /// \return true = success; false = failure.
//    bool OpenForWrite(const char * filename);
//
//    /// Reads the specified bam index file.  It must be read prior to setting a
//    /// read section, for seeking and reading portions of a bam file.
//    /// \return true = success; false = failure.   
//    bool ReadBamIndex(const char * filename);
//
//    /// Close the file if there is one open.
//    void Close();
//
//    /// Returns whether or not the end of the file has been reached.
//    /// \return true = EOF; false = not eof.
//    /// If the file is not open, false is returned.
//    bool IsEOF();
//   
//    /// Reads the header section from the file and stores it in
//    /// the passed in header.
//    /// \return true = success; false = failure.
//    bool ReadHeader(SamFileHeader& header);
//   
//    /// Writes the specified header into the file.
//    /// \return true = success; false = failure.
//    bool WriteHeader(SamFileHeader& header);
//
//    /// Reads the next record from the file & stores it in the passed in record.
//    /// \return true  = record was successfully set.
//    ///                false = record was not successfully set.
//    bool ReadRecord(SamFileHeader& header, SamRecord& record);
//   
//    /// Writes the specified record into the file.
//    /// \return true = success; false = failure.
//    bool WriteRecord(SamFileHeader& header, SamRecord& record);
//   
//    /// Set the flag to validate that the file is sorted as it is read/written.
//    /// Must be called after the file has been opened.
//    void setSortedValidation(SortedType sortType);
//
//    /// Return the number of records that have been read/written so far.
//    uint32_t GetCurrentRecordCount();
//
//    /// Get the Status of the last call that sets status.
//    /// To remain backwards compatable - will be removed later.
//    inline SamStatus::Status GetFailure()
//    {
//        return(GetStatus());
//    }
//
//    /// Get the Status of the last call that sets status.
//    inline SamStatus::Status GetStatus()
//    {
//        return(myStatus.getStatus());
//    }
//
//    /// Get the Status of the last call that sets status.
//    inline const char* GetStatusMessage()
//    {
//        return(myStatus.getStatusMessage());
//    }
//
//    /// Sets what part of the BAM file should be read.  This version will
//    /// set it to only read a specific reference id.  The records for that
//    /// reference id will be retrieved on each ReadRecord call.  When all
//    /// records have been retrieved for the specified reference id, ReadRecord
//    /// will return failure until a new read section is set.
//    /// Must be called only after the file has been opened for reading.
//    /// \param  refID the reference ID of the records to read from the file.
//    /// \return true = success; false = failure.
//    bool SetReadSection(int32_t refID);
//
//    /// Sets what part of the BAM file should be read.  This version will
//    /// set it to only read a specific reference name.  The records for that
//    /// reference id will be retrieved on each ReadRecord call.  When all
//    /// records have been retrieved for the specified reference name,
//    /// ReadRecord will return failure until a new read section is set.
//    /// Must be called only after the file has been opened for reading.
//    /// \param  refName the reference name of the records to read from the file.
//    /// \return true = success; false = failure.
//    bool SetReadSection(const char* refName);
//
//    /// Sets what part of the BAM file should be read.  This version will
//    /// set it to only read a specific reference id and start/end position.
//    /// The records for this section will be retrieved on each ReadRecord
//    /// call.  When all records have been retrieved for the specified section,
//    /// ReadRecord will return failure until a new read section is set.
//    /// Must be called only after the file has been opened for reading.
//    /// \param  refID the reference ID of the records to read from the file.
//    /// \param  start inclusive start position of records that should be read for this refID.
//    /// \param  end exclusive end position of records that should be read for this refID.
//    /// \return true = success; false = failure.   
//    bool SetReadSection(int32_t refID, int32_t start, int32_t end);
//
//    /// Sets what part of the BAM file should be read.  This version will
//    /// set it to only read a specific reference name and start/end position.
//    /// The records for this section will be retrieved on each ReadRecord
//    /// call.  When all records have been retrieved for the specified section,
//    /// ReadRecord will return failure until a new read section is set.
//    /// Must be called only after the file has been opened for reading.
//    /// \param  refName the reference name of the records to read from the file.
//    /// \param  start inclusive start position of records that should be read for this refID.
//    /// \param  end exclusive end position of records that should be read for this refID.
//    /// \return true = success; false = failure.   
//    bool SetReadSection(const char* refName, int32_t start, int32_t end);
//
//    /// Whether or not statistics should be generated for this file.
//    /// The value is carried over between files and is not reset, but
//    /// the statistics themselves are reset between files.
//    /// \param genStats set to true if statistics should be generated, false if not.
//    void GenerateStatistics(bool genStats);
//
//    inline void PrintStatistics() {if(myStatistics != NULL) myStatistics->print();}
//
//protected:
//    void resetFile();
//
//    /// Validate that the record is sorted compared to the previously read record
//    /// if there is one, according to the specified sort order.
//    /// If the sort order is UNSORTED, true is returned.
//    bool validateSortOrder(SamRecord& record, SamFileHeader& header);
//   
//    // Return the sort order as defined by the header.  If it is undefined
//    // or set to an unknown value, UNSORTED is returned.
//    SortedType getSortOrderFromHeader(SamFileHeader& header);
//
//    /// Overwrites read record to read from the specific reference only.
//    bool readIndexedRecord(SamFileHeader& header, SamRecord& record);
//
//    bool processNewSection(SamFileHeader &header);
//
//    IFILE  myFilePtr;
//    GenericSamInterface* myInterfacePtr;
//
//    /// Flag to indicate if a file is open for reading.
//    bool myIsOpenForRead;
//    /// Flag to indicate if a file is open for writing.
//    bool myIsOpenForWrite;
//    /// Flag to indicate if a header has been read/written - required before
//    /// being able to read/write a record.
//    bool myHasHeader;
//
//    SortedType mySortedType;
//
//    /// Previous values used for checking if the file is sorted.
//    int32_t myPrevCoord;
//    int32_t myPrevRefID;
//    std::string myPrevReadName;
//
//    /// Keep a count of the number of records that have been read/written so far.
//    uint32_t myRecordCount;
//
//    /// Pointer to the statistics for this file.
//    SamStatistics* myStatistics;
//   
//    /// The status of the last SamFile command.
//    SamStatus myStatus;
//
//    /// Values for reading Sorted BAM files via the index.
//    bool myIsBamOpenForRead;
//    bool myNewSection;
//    int32_t myRefID;
//    int32_t myStartPos;
//    int32_t myEndPos;
//    uint64_t myCurrentChunkEnd;
//    SortedChunkList myChunksToRead;
//    BamIndex* myBamIndex;
//
//    std::string myRefName;
//};
//
//
//class SamFileReader : public SamFile
//{
//public:
//
//    /// Default Constructor.
//    SamFileReader();
//
//    /// Constructor that opens the specified file for read.
//    SamFileReader(const char* filename);
//
//    virtual ~SamFileReader();
//};
//
//
//class SamFileWriter : public SamFile
//{
//public:
//    /// Default Constructor.
//    SamFileWriter();
//
//    /// Constructor that opens the specified file for write.
//    SamFileWriter(const char* filename);
//
//    virtual ~SamFileWriter();
//};
//
//#endif
//
//
//
//
//#ifndef __SAM_VALIDATION_H__
//#define __SAM_VALIDATION_H__
//
//#include "SamFile.h"
//#include <list>
//
//class SamValidationError
//{
//public:
//    // Warning is used if it is just an invalid value.
//    // Error is used if parsing could not succeed.
//    enum Severity {WARNING, ERROR};
//
//    enum Type
//        {
//            // 
//            INVALID_QNAME,
//            INVALID_REF_ID,
//            INVALID_RNAME,
//            INVALID_POS,
//            INVALID_MAPQ,
//            INVALID_CIGAR,
//            INVALID_MRNM,
//            INVALID_QUAL
//        };
//
//    static const char* getTypeString(Type type);
//
//    SamValidationError(Type type, Severity severity, std::string Message);
//   
//    Type getType() const;
//    Severity getSeverity() const;
//    const char* getMessage() const;
//
//    const char* getTypeString() const;
//    const char* getSeverityString() const;
//
//    void getErrorString(std::string& errorString) const;
//
//    void printError() const;
//
//private:
//    SamValidationError();
//
//    static const char* enumTypeString[];
//    static const char* enumSeverityString[];
//
//    Type myType;
//    Severity mySeverity;
//    std::string myMessage;
//
//};
//
//
////
//// stream output for validation failure information
////
//inline std::ostream &operator << (std::ostream &stream, 
//                                  const SamValidationError &error)
//{
//    std::string errorMessage;
//    error.getErrorString(errorMessage);
//    stream << errorMessage;
//    return stream;
//}
//
//
//class SamValidationErrors
//{
//public:
//    // Constructor.
//    SamValidationErrors();
//    // Destructor
//    ~SamValidationErrors();
//
//    // Remove all the errors from the list.
//    void clear();
//
//    // Adds the specified error to the set of errors.
//    void addError(SamValidationError::Type newType, 
//                  SamValidationError::Severity newSeverity,
//                  const char* newMessage);
//
//    // Return the number of validation errors that are contained in this object.
//    unsigned int numErrors();
//
//    // Return a pointer to the next error.  It does not remove it from the list.
//    // Returns null once all errors have been retrieved until resetErrorIter
//    // is called.
//    const SamValidationError* getNextError();
//   
//    // Resets the iterator to the begining of the errors.
//    void resetErrorIter();
//
//    // Appends the error messages to the passed in string.
//    void getErrorString(std::string& errorString) const;
//
//private:
//    std::list<const SamValidationError*> myValidationErrors;
//    std::list<const SamValidationError*>::const_iterator myErrorIter;
//};
//
//
////
//// stream output for all validation failures information
////
//inline std::ostream& operator << (std::ostream& stream,
//                                  const SamValidationErrors& errors)
//{
//    std::string errorString = "";
//    errors.getErrorString(errorString);
//    stream << errorString;
//    return stream;
//}
//
//
//class SamValidator
//{
//public:
//
//    static bool isValid(SamFileHeader& samHeader, SamRecord& samRecord, 
//                        SamValidationErrors& validationErrors);
//
//    static bool isValidQname(const char* qname, uint8_t qnameLen, 
//                             SamValidationErrors& validationErrors);
//    static bool isValidFlag(uint16_t flag,
//                            SamValidationErrors& validationErrors);
//    // Validate the rname including validating against the header.
//    static bool isValidRname(SamFileHeader& samHeader, 
//                             const char* rname,
//                             SamValidationErrors& validationErrors);
//    // Validate the rname without validating against the header.
//    static bool isValidRname(const char* rname,
//                             SamValidationErrors& validationErrors);
//    static bool isValidRefID(int32_t refID, const StringArray& refContigs, 
//                             SamValidationErrors& validationErrors);
//    static bool isValid1BasedPos(int32_t pos, 
//                                 SamValidationErrors& validationErrors);
//    static bool isValidMapQuality(uint8_t mapQuality,
//                                  SamValidationErrors& validationErrors);
//    // Cigar validation depends on sequence.
//    static bool isValidCigar(const char* cigar, const char* sequence,
//                             SamValidationErrors& validationErrors);
//    static bool isValidMrnm();
//    static bool isValidMpos();
//    static bool isValidIsize();
//    static bool isValidSeq();
//    // Quality validation depends on sequence.
//    static bool isValidQuality(const char* quality, const char* sequence,
//                                  SamValidationErrors& validationErrors);
//    static bool isValidTag();
//    static bool isValidVtype();
//    static bool isValidValue();
//};
//
//
//#endif


#include "SamFile.h"
#include "SamValidation.h"



int ReadIndexedBam(const char* inputFilename,
                   const char* outputFilename,
                   const char* indexFilename)
{

   //// Open the input file for reading.
   SamFile samIn;
   //samIn.OpenForRead(inputFilename);
   //
   //// Open the bam index file for reading.
   //samIn.ReadBamIndex(indexFilename);
   //
   //// Open the output file for writing.
   //SamFile samOut;
   //samOut.OpenForWrite(outputFilename);
   //
   //// Read the sam header.
   //SamFileHeader samHeader;
   //samIn.ReadHeader(samHeader);
   //
   //// Write the sam header.
   //samOut.WriteHeader(samHeader);
   //
   //SamRecord samRecord;
   //SamValidationErrors samInvalidErrors;
   //
   //// Loop through each Reference.
   //for(int i = -1; i < 23; i++)
   //{
   //   int numSectionRecords = 0;
   //   samIn.SetReadSection(i);
   //
   //   // Keep reading records until they aren't more.
   //   while(samIn.ReadRecord(samHeader, samRecord))
   //   {
   //      numSectionRecords++;
   //
   //      // Successfully read a record from the file, so write it.
   //      samOut.WriteRecord(samHeader, samRecord);
   //   }
   //
   //   std::cout << "Reference ID " << i << " has " << numSectionRecords 
   //             << " records" << std::endl;
   //}
   //
   //std::cout << "Number of records = " << samIn.GetCurrentRecordCount() << std::endl;

   return(0);
}


int main(int argc, char *argv[]) {

   const char* inputFilename;
   const char* outputFilename;
   const char* indexFilename;

   ReadIndexedBam(inputFilename, outputFilename, indexFilename);
}


