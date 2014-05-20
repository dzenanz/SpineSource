#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <string>
#include <iterator>
#include "itkDatImageIO.h"
#include "itkMacro.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include "itksys/SystemTools.hxx"
#include <fstream>
#include <iomanip>
#include <ios>

typedef std::string::value_type char_t;

char_t up_char( char_t ch )
{
    return std::use_facet< std::ctype< char_t > >( std::locale() ).toupper( ch );
}

std::string toupper( const std::string &src )
{
    std::string result;
    std::transform( src.begin(), src.end(), std::back_inserter( result ), up_char );
    return result;
}

namespace itk {

	bool DatImageIO::SupportsDimension(unsigned long dim)
	{
		return dim==3;
	}

	void DatImageIO::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	bool DatImageIO::CanReadFile( const char* filename ) 
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
		std::string::size_type datPos = fname.rfind(".dat");
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

	void DatImageIO::ReadImageInformation()
	{
		std::ifstream in(this->GetFileName());

		std::string s, localRawFilename;
		char buf[1000];
		in>>s; in>>localRawFilename;
		in.getline(buf, 999); //read endline char
		in.getline(buf, 999); //read TaggedFileName line

		int sepPos=localRawFilename.find_last_of("/\\");
		if (sepPos==localRawFilename.npos)
		{
			s=this->GetFileName();
			sepPos=s.find_last_of("/\\");
			if (sepPos!=s.npos)
			{
				s.resize(sepPos+1);
				rawFilename=s+localRawFilename;
			}
			else
				rawFilename=localRawFilename;
		}
		else
			rawFilename=localRawFilename;

		this->SetNumberOfDimensions(3);
		this->SetPixelType( ImageIOBase::SCALAR );
		this->SetNumberOfComponents(1);

		int xres, yres, zres;
		in>>s; in>>xres>>yres>>zres;
		in.getline(buf, 999); //read endline char
		this->SetDimensions(0, xres);
		this->SetDimensions(1, yres);
		this->SetDimensions(2, zres);

		double xs, ys, zs;
		in>>s; in>>xs>>ys>>zs;
		in.getline(buf, 999); //read endline char
		this->SetSpacing(0, xs);
		this->SetSpacing(1, ys);
		this->SetSpacing(2, zs);

		in>>s; in>>s;
		in.close();

		if (s=="UCHAR")
			this->SetComponentType(UCHAR);
		else if (s=="CHAR")
			this->SetComponentType(CHAR);
		else if (s=="USHORT")
			this->SetComponentType(USHORT);
		else if (s=="SHORT")
			this->SetComponentType(SHORT);
		else if (s=="USHORT_12BIT")
			{
				ExceptionObject exception(__FILE__, __LINE__);
				exception.SetDescription("12-bit data is not supported currently.");
				throw exception;
			}
		else if (s=="UINT")
			this->SetComponentType(UINT);
		else if (s=="INT")
			this->SetComponentType(INT);
		else if (s=="ULONG")
			this->SetComponentType(ULONG);
		else if (s=="LONG")
			this->SetComponentType(LONG);
		else if (s=="FLOAT")
			this->SetComponentType(FLOAT);
		else if (s=="DOUBLE")
			this->SetComponentType(DOUBLE);
		else
			{
				ExceptionObject exception(__FILE__, __LINE__);
				exception.SetDescription("Format type not recognized.");
				throw exception;
			}
	}


	void DatImageIO::Read(void* buffer)
	{
		size_t elems=this->GetDimensions(0)*this->GetDimensions(1)*this->GetDimensions(2);
		size_t esize=this->GetComponentSize();
		FILE *f=fopen(rawFilename.c_str(), "rb");
		if (fread(buffer, esize, elems, f)!=elems)
		{
			ExceptionObject exception(__FILE__, __LINE__);
			exception.SetDescription("Error executing fread on "+rawFilename);
			throw exception;
		}
		fclose(f);		
	} 


	bool DatImageIO::CanWriteFile( const char * name )
	{
		std::string filename = name;
		if(  filename == "" )
		{
			return false;
		}

		std::string::size_type datPos = filename.rfind(".dat");
		if ((datPos != std::string::npos)
			&& (datPos == filename.length() - 4))
		{
			return true;
		}

		return false;
	}

	void DatImageIO::WriteImageInformation(void)
	{
        std::string filename=this->GetFileName();
        std::ofstream f(this->GetFileName());
        f.setf(std::ios::left);
        if (!f.good())
		{
			ExceptionObject exception(__FILE__, __LINE__);
            exception.SetDescription("Could not open file for writing: "+filename);
			throw exception;
		}
        filename=itksys::SystemTools::GetFilenameWithoutLastExtension(filename)+".raw";
        f<<std::setw(16)<<"ObjectFileName:"<<filename<<std::endl;
        f<<"TaggedFileName: ---\n";
        f<<std::setw(16)<<"Resolution:"<<this->GetDimensions(0)<<' '<<this->GetDimensions(1)<<' '<<this->GetDimensions(2)<<std::endl;
        f<<std::setw(16)<<"SliceThickness:"<<this->GetSpacing(0)<<' '<<this->GetSpacing(1)<<' '<<this->GetSpacing(2)<<std::endl;
        f<<std::setw(16)<<"Format:"<<toupper(this->GetComponentTypeAsString(this->GetComponentType()))<<std::endl;
        f<<"NbrTags:        0\nObjectType:     TEXTURE_VOLUME_OBJECT\nObjectModel:    RGBA\nGridType:       EQUIDISTANT";
        f.close();
    }


	void DatImageIO::Write( const void* buffer) 
	{
        WriteImageInformation();
        std::string filename=itksys::SystemTools::GetFilenamePath(this->GetFileName())+"/"
            +itksys::SystemTools::GetFilenameWithoutLastExtension(this->GetFileName())+".raw";
        size_t elems=this->GetDimensions(0)*this->GetDimensions(1)*this->GetDimensions(2);
		FILE *f=fopen(filename.c_str(), "wb");
		if (fwrite(buffer, this->GetComponentSize(), elems, f)!=elems)
		{
			ExceptionObject exception(__FILE__, __LINE__);
			exception.SetDescription("Error executing fwrite on "+filename);
			throw exception;
		}
		fclose(f);		
	}

} // end namespace itk
