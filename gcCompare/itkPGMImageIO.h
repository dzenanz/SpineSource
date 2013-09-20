#ifndef __itkPGMImageIO_h
#define __itkPGMImageIO_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkImageIOBase.h"
#include <fstream>

namespace itk
{

	/** \class PGMImageIO
	*
	* \brief Read and write the Netpbm format's PGM extension to 3D.
	*
	*  \ingroup IOFilters
	*/
	class ITK_EXPORT PGMImageIO : public ImageIOBase
	{
	public:
		/** Standard class typedefs. */
		typedef PGMImageIO        Self;
		typedef ImageIOBase        Superclass;
		typedef SmartPointer<Self> Pointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(PGMImageIO, ImageIOBase);

		/** The different types of ImageIO's can support data of varying
		* dimensionality. For example, some file formats are strictly 2D
		* while others can support 2D, 3D, or even n-D. This method returns
		* true/false as to whether the ImageIO can support the dimension
		* indicated. */
		virtual bool SupportsDimension(unsigned long );

		/** Determine the file type. Returns true if this ImageIO can read the
		* file specified. */
		virtual bool CanReadFile(const char*);

		/** Set the spacing and dimension information for the set filename. */
		virtual void ReadImageInformation();

		/** Reads the data from disk into the memory buffer provided. */
		virtual void Read(void* buffer);

		/** Determine the file type. Returns true if this ImageIO can write the
		* file specified. */
		virtual bool CanWriteFile(const char*);

		virtual void WriteImageInformation();

		/** Writes the data to disk from the memory buffer provided. Make sure
		* that the IORegions has been set properly. */
		virtual void Write(const void* buffer);

	protected:
		PGMImageIO() {}
		~PGMImageIO() {}
		void PrintSelf(std::ostream& os, Indent indent) const;

	private:
		PGMImageIO(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented

	};

} // end namespace itk

#endif // __itkPGMImageIO_h
