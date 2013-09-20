#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <string>
#include <iterator>
#include "itkPGMImageIO.h"
#include "itkMacro.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "itksys/SystemTools.hxx"
#include <fstream>
#include <iomanip>
#include <ios>

typedef std::string::value_type char_t;

namespace itk {

	bool PGMImageIO::SupportsDimension(unsigned long dim)
	{
		return dim==3 || dim==2; //2D + extension to 3D
	}

	void PGMImageIO::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	bool PGMImageIO::CanReadFile( const char* filename ) 
	{
		// Check the extension first to avoid opening files that do not
		// look like dats.  The file must have an appropriate extension to be
		// recognized.
		std::string fname = filename;

		if(  fname == "" )
		{
			itkDebugMacro(<<"No filename specified.");
			return false;
		}

		bool extensionFound = false;
		std::string::size_type datPos = fname.rfind(".pgm");
		if ((datPos != std::string::npos)
			&& (datPos == fname.length() - 4))
		{
			extensionFound = true;
		}

		if( !extensionFound )
		{
			itkDebugMacro(<<"The filename extension is not recognized");
			return false;
		}

		std::ifstream inputStream;

		inputStream.open( filename, std::ios::in | std::ios::binary );

		if( inputStream.fail() )
		{
			return false;
		}

		inputStream.close();
		return true;
	}

	void PGMImageIO::ReadImageInformation()
	{
		std::ifstream in(this->GetFileName());

		std::string s;
		char buf[1000];
		in>>s;
		if (s!="P5")
		{
			ExceptionObject exception(__FILE__, __LINE__);
			exception.SetDescription("P5 signature not found.");
			throw exception;
		}
		in.getline(buf, 999); //read endline char

		int xres, yres, zres, maxVal;
		in>>xres>>yres;

        this->SetComponentType(UCHAR);
		this->SetPixelType( ImageIOBase::SCALAR );
		this->SetNumberOfComponents(1);

        in.getline(buf, 999);
        if (buf[0]!=0)
        {
            this->SetNumberOfDimensions(3);
            zres=atoi(buf);
		    this->SetDimensions(2, zres);
        }
        else
        {
            this->SetNumberOfDimensions(2);
        }
		this->SetDimensions(0, xres);
		this->SetDimensions(1, yres);

        in>>maxVal;
        if (maxVal!=255)
            itkWarningMacro (<<"Maximum value should be 255, not "<<maxVal);
		in.close();
	}

	void PGMImageIO::Read(void* buffer)
	{
        size_t elems;
        if (this->GetNumberOfDimensions()==3)
            elems=this->GetDimensions(0)*this->GetDimensions(1)*this->GetDimensions(2);
        else
            elems=this->GetDimensions(0)*this->GetDimensions(1);
		size_t esize=this->GetComponentSize();
        std::ifstream in(this->GetFileName(), std::ios::binary);
        in.getline((char *)buffer, 999);
        in.getline((char *)buffer, 999);
        in.getline((char *)buffer, 999);
        in.read((char *)buffer, esize*elems);
        if (in.bad()||in.gcount()!=esize*elems)
		{
			ExceptionObject exception(__FILE__, __LINE__);
			exception.SetDescription("Error reading "+std::string(this->GetFileName()));
			throw exception;
		}
		in.close();
	} 

	bool PGMImageIO::CanWriteFile( const char * name )
	{
		std::string filename = name;
		if(  filename == "" )
		{
			return false;
		}

		std::string::size_type datPos = filename.rfind(".pgm");
		if ((datPos != std::string::npos)
			&& (datPos == filename.length() - 4))
            return true;

		return false;
	}

	void PGMImageIO::WriteImageInformation(void)
	{
        //nothing to do here
    }

	void PGMImageIO::Write( const void* buffer) 
	{
        size_t elems;
        if (this->GetNumberOfDimensions()==3)
            elems=this->GetDimensions(0)*this->GetDimensions(1)*this->GetDimensions(2);
        else
            elems=this->GetDimensions(0)*this->GetDimensions(1);
		std::ofstream f(this->GetFileName(), std::ios::binary);
        f<<"P5\n";
        f<<this->GetDimensions(0)<<' '<<this->GetDimensions(1);
        if (this->GetNumberOfDimensions()==3)
            f<<' '<<this->GetDimensions(2);
        f<<"\n255\n";
        f.write((char *)buffer, elems*this->GetComponentSize());

        if (f.bad())
		{
			ExceptionObject exception(__FILE__, __LINE__);
			exception.SetDescription("Error writing "+std::string(this->GetFileName()));
			throw exception;
		}
		f.close();		
	}

} // end namespace itk
