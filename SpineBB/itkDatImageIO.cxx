#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <string>
#include "itkDatImageIO.h"
#include "itkMacro.h"
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"
#include <fstream>

#if defined(__BORLANDC__) 
# include <math.h> 
# include <float.h> // for _control87() 
#endif // defined(__BORLANDC__) 

namespace itk {

	bool DatImageIO::SupportsDimension(unsigned long dim)
	{
		return dim==3;
	}

	void DatImageIO::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
	}

	//ImageIOBase::IOComponentType
	//DatImageIO::
	//DatToITKComponentType( const int datComponentType ) const
	//{
	//#if defined(__BORLANDC__) 
	//// Disable floating point exceptions in Borland 
	//  _control87(MCW_EM, MCW_EM); 
	//#endif // defined(__BORLANDC__) 
	//  switch( datComponentType )
	//    {
	//    case datTypeUnknown:
	//    case datTypeBlock:
	//      return UNKNOWNCOMPONENTTYPE;
	//      break;
	//    case datTypeChar:
	//      return CHAR;
	//      break;
	//    case datTypeUChar:
	//      return UCHAR;
	//      break;
	//    case datTypeShort:
	//      return SHORT;
	//      break;
	//    case datTypeUShort:
	//      return USHORT;
	//      break;
	//    // "long" is a silly type because it basically guaranteed not to be
	//    // cross-platform across 32-vs-64 bit machines, but we'll use it 
	//    // where possible.
	//    case datTypeLLong:
	//      return airMy32Bit ? UNKNOWNCOMPONENTTYPE : LONG;
	//      break;
	//    case datTypeULLong:
	//      return airMy32Bit ? UNKNOWNCOMPONENTTYPE : ULONG;
	//      break;
	//    case datTypeInt:
	//      return INT;
	//      break;
	//    case datTypeUInt:
	//      return UINT;
	//      break;
	//    case datTypeFloat:
	//      return FLOAT;
	//      break;
	//    case datTypeDouble:
	//      return DOUBLE;
	//      break;
	//    default:
	//      return UNKNOWNCOMPONENTTYPE;
	//      break;
	//    }
	//
	//  // Strictly to avoid compiler warning regarding "control may reach end of
	//  // non-void function":
	//  //
	//  return UNKNOWNCOMPONENTTYPE;
	//}
	//
	//int
	//DatImageIO::
	//ITKToDatComponentType( const ImageIOBase::IOComponentType itkComponentType ) const
	//{
	//#if defined(__BORLANDC__) 
	//// Disable floating point exceptions in Borland 
	//  _control87(MCW_EM, MCW_EM); 
	//#endif // defined(__BORLANDC__) 
	//  switch( itkComponentType )
	//    {
	//    case UNKNOWNCOMPONENTTYPE:
	//      return datTypeUnknown;
	//      break;
	//    case CHAR:
	//      return datTypeChar;
	//      break;
	//    case UCHAR:
	//      return datTypeUChar;
	//      break;
	//    case SHORT:
	//      return datTypeShort;
	//      break;
	//    case USHORT:
	//      return datTypeUShort;
	//      break;
	//    // "long" is a silly type because it basically guaranteed not to be
	//    // cross-platform across 32-vs-64 bit machines, but we can figure out
	//    // a cross-platform way of storing the information.
	//    case LONG:
	//      return airMy32Bit ? datTypeInt : datTypeLLong;
	//      break;
	//    case ULONG:
	//      return airMy32Bit ? datTypeUInt : datTypeULLong;
	//      break;
	//    case INT:
	//      return datTypeInt;
	//      break;
	//    case UINT:
	//      return datTypeUInt;
	//      break;
	//    case FLOAT:
	//      return datTypeFloat;
	//      break;
	//    case DOUBLE:
	//      return datTypeDouble;
	//      break;
	//    default:
	//      return datTypeUnknown;
	//      break;
	//    }
	//
	//  // Strictly to avoid compiler warning regarding "control may reach end of
	//  // non-void function":
	//  //
	//  return datTypeUnknown;
	//}

	bool DatImageIO::CanReadFile( const char* filename ) 
	{
#if defined(__BORLANDC__) 
		// Disable floating point exceptions in Borland 
		_control87(MCW_EM, MCW_EM); 
#endif // defined(__BORLANDC__) 
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
#if defined(__BORLANDC__) 
		// Disable floating point exceptions in Borland 
		_control87(MCW_EM, MCW_EM); 
#endif // defined(__BORLANDC__) 

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
		else if (s=="USHORT")
			this->SetComponentType(USHORT);
		else if (s=="USHORT_12BIT")
			{
				ExceptionObject exception(__FILE__, __LINE__);
				exception.SetDescription("12-bit data is not supported currently.");
				throw exception;
			}
		else if (s=="UINT")
			this->SetComponentType(UINT);
		else if (s=="ULONG")
			this->SetComponentType(ULONG);
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
#if defined(__BORLANDC__) 
		// Disable floating point exceptions in Borland 
		_control87(MCW_EM, MCW_EM); 
#endif // defined(__BORLANDC__) 

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
#if defined(__BORLANDC__) 
		// Disable floating point exceptions in Borland 
		_control87(MCW_EM, MCW_EM); 
#endif // defined(__BORLANDC__) 

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
		// Nothing needs doing here.
	}


	void DatImageIO::Write( const void* buffer) 
	{
		throw itk::ExceptionObject(__FILE__, __LINE__, "Writing Erlangen format (.dat) is not implemented", ".dat export filter");

		//to implement

		//Dat *dat = datNew();
		//DatIoState *nio = datIoStateNew();
		//int kind[DAT_DIM_MAX];
		//size_t size[DAT_DIM_MAX];
		//unsigned int datDim, baseDim, spaceDim;
		//double spaceDir[DAT_DIM_MAX][DAT_SPACE_DIM_MAX];
		//double origin[DAT_DIM_MAX];

		//spaceDim = this->GetNumberOfDimensions();
		//if (this->GetNumberOfComponents() > 1)
		//  {
		//  size[0] = this->GetNumberOfComponents();
		//  switch (this->GetPixelType())
		//    {
		//    case ImageIOBase::RGB:
		//      kind[0] = datKindRGBColor;
		//      break;
		//    case ImageIOBase::RGBA:
		//      kind[0] = datKindRGBAColor;
		//      break;
		//    case ImageIOBase::POINT:
		//      kind[0] = datKindPoint;
		//      break;
		//    case ImageIOBase::COVARIANTVECTOR:
		//      kind[0] = datKindCovariantVector;
		//      break;
		//    case ImageIOBase::SYMMETRICSECONDRANKTENSOR:
		//    case ImageIOBase::DIFFUSIONTENSOR3D:
		//      kind[0] = datKind3DSymMatrix;
		//      break;
		//    case ImageIOBase::COMPLEX:
		//      kind[0] = datKindComplex;
		//      break;
		//    case ImageIOBase::VECTOR:
		//    case ImageIOBase::OFFSET:      // HEY is this right?
		//    case ImageIOBase::FIXEDARRAY:  // HEY is this right?
		//    default:
		//      kind[0] = datKindVector;
		//      break;
		//    }
		//  // the range axis has no space direction
		//  for (unsigned int saxi=0; saxi < spaceDim; saxi++)
		//    {
		//    spaceDir[0][saxi] = AIR_NAN;
		//    }
		//  baseDim = 1;
		//  }
		//else
		//  {
		//  baseDim = 0;
		//  }
		//datDim = baseDim + spaceDim;
		//std::vector<double> spaceDirStd(spaceDim);
		//unsigned int axi;
		//for (axi=0; axi < spaceDim; axi++)
		//  {
		//  size[axi+baseDim] = this->GetDimensions(axi);
		//  kind[axi+baseDim] = datKindDomain;
		//  origin[axi] = this->GetOrigin(axi);
		//  double spacing = this->GetSpacing(axi);
		//  spaceDirStd = this->GetDirection(axi);
		//  for (unsigned int saxi=0; saxi < spaceDim; saxi++)
		//    {
		//    spaceDir[axi+baseDim][saxi] = spacing*spaceDirStd[saxi];
		//    }
		//  }
		//if (datWrap_nva(dat, const_cast<void *>(buffer),
		//                 this->ITKToDatComponentType( m_ComponentType ),
		//                 datDim, size) || (3 == spaceDim
		//     // special case: ITK is LPS in 3-D
		//     ? datSpaceSet(dat, datSpaceLeftPosteriorSuperior)
		//     : datSpaceDimensionSet(dat, spaceDim)) ||
		//    datSpaceOriginSet(dat, origin))
		//  {
		//  char *err = biffGetDone(DAT); // would be nice to free(err)
		//  itkExceptionMacro("Write: Error wrapping dat for " 
		//                    << this->GetFileName() << ":\n" << err);
		//  }
		//datAxisInfoSet_nva(dat, datAxisInfoKind, kind);
		//datAxisInfoSet_nva(dat, datAxisInfoSpaceDirection, spaceDir);

		//// Go through MetaDataDictionary and set either specific dat field
		//// or a key/value pair
		//MetaDataDictionary &thisDic = this->GetMetaDataDictionary();
		//std::vector<std::string> keys = thisDic.GetKeys();
		//std::vector<std::string>::const_iterator keyIt;
		//const char *keyField, *field;
		//for( keyIt = keys.begin(); keyIt != keys.end(); keyIt++ )
		//  {
		//  if (!strncmp(KEY_PREFIX, (*keyIt).c_str(), strlen(KEY_PREFIX)))
		//    {
		//    keyField = (*keyIt).c_str() + strlen(KEY_PREFIX);
		//    // only of one of these can succeed
		//    field = airEnumStr(datField, datField_thicknesses);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      if (1 == sscanf(keyField + strlen(field), "[%d]", &axi)
		//          && axi + baseDim < dat->dim)
		//        {
		//        double thickness = 0.0;  // local for Borland
		//        ExposeMetaData<double>(thisDic, *keyIt, thickness);
		//        dat->axis[axi+baseDim].thickness = thickness;
		//        }
		//      }
		//    field = airEnumStr(datField, datField_centers);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      if (1 == sscanf(keyField + strlen(field), "[%d]", &axi)
		//          && axi + baseDim < dat->dim)
		//        {
		//        std::string value;  // local for Borland
		//        ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//        dat->axis[axi+baseDim].center = airEnumVal(datCenter,
		//                                                    value.c_str());
		//        }
		//      }
		//    field = airEnumStr(datField, datField_kinds);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      if (1 == sscanf(keyField + strlen(field), "[%d]", &axi)
		//          && axi + baseDim < dat->dim)
		//        {
		//        std::string value;  // local for Borland
		//        ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//        dat->axis[axi+baseDim].kind = airEnumVal(datKind,
		//                                                  value.c_str());
		//        }
		//      }
		//    field = airEnumStr(datField, datField_labels);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      if (1 == sscanf(keyField + strlen(field), "[%d]", &axi)
		//          && axi + baseDim < dat->dim)
		//        {
		//        std::string value;  // local for Borland
		//        ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//        dat->axis[axi+baseDim].label = airStrdup(value.c_str());
		//        }
		//      }
		//    field = airEnumStr(datField, datField_old_min);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      ExposeMetaData<double>(thisDic, *keyIt, dat->oldMin);
		//      }
		//    field = airEnumStr(datField, datField_old_max);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      ExposeMetaData<double>(thisDic, *keyIt, dat->oldMax);
		//      }

		//    field = airEnumStr(datField, datField_space);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      int space;
		//      std::string value;  // local for Borland
		//      ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//      space = airEnumVal(datSpace, value.c_str());
		//      if (datSpaceDimension(space) == dat->spaceDim)
		//        {
		//        // sanity check
		//        dat->space = space;
		//        }
		//      }

		//    field = airEnumStr(datField, datField_content);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      std::string value;  // local for Borland
		//      ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//      dat->content = airStrdup(value.c_str());
		//      }
		//    field = airEnumStr(datField, datField_measurement_frame);
		//    if (!strncmp(keyField, field, strlen(field)))
		//      {
		//      std::vector<std::vector<double> > msrFrame;
		//      ExposeMetaData<std::vector<std::vector<double> > >(thisDic,
		//                                                         *keyIt, msrFrame);
		//      for (unsigned int saxi=0; saxi < dat->spaceDim; saxi++)
		//        {
		//        for (unsigned int saxj=0; saxj < dat->spaceDim; saxj++)
		//          {
		//          if (saxi < msrFrame.size() &&
		//              saxj < msrFrame[saxi].size())
		//            {
		//            dat->measurementFrame[saxi][saxj] = msrFrame[saxi][saxj];
		//            }
		//          else
		//            {
		//            // there is a difference between the dimension of the 
		//            // recorded measurement frame, and the actual dimension of
		//            // the ITK image, which (for now) determines dat->spaceDim.
		//            // We can't set this to AIR_NAN, because the coefficients of
		//            // the measurement frame have to all be equally existent.
		//            // If we used 0, it might not a flag that something is wrong.
		//            // So, we have to get creative.
		//            dat->measurementFrame[saxi][saxj] = 666666;
		//            }
		//          }
		//        }
		//      }
		//    }
		//  else
		//    {
		//    // not a DAT field packed into meta data; just a regular key/value
		//    std::string value;  // local for Borland
		//    ExposeMetaData<std::string>(thisDic, *keyIt, value);
		//    datKeyValueAdd(dat, (*keyIt).c_str(), value.c_str());
		//    }
		//  }

		//// set encoding for data: compressed (raw), (uncompressed) raw, or ascii
		//if (this->GetUseCompression() == true
		//    && datEncodingGzip->available())
		//  {
		//  // this is necessarily gzip-compressed *raw* data
		//  nio->encoding = datEncodingGzip;
		//  }
		//else
		//  {
		//  Superclass::FileType fileType = this->GetFileType();
		//  switch ( fileType )
		//    {
		//    default:
		//    case TypeNotApplicable:
		//    case Binary:
		//      nio->encoding = datEncodingRaw;
		//      break;
		//    case ASCII:
		//      nio->encoding = datEncodingAscii;
		//      break;
		//    }
		//  }

		//// set desired endianness of output
		//Superclass::ByteOrder byteOrder = this->GetByteOrder();
		//switch (byteOrder)
		//  {
		//  default:
		//  case OrderNotApplicable:
		//    nio->endian = airEndianUnknown;
		//    break;
		//  case BigEndian:
		//    nio->endian = airEndianBig;
		//    break;
		//  case LittleEndian:
		//    nio->endian = airEndianLittle;
		//    break;
		//  }
		//
		//// Write the dat to file.
		//if (datSave(this->GetFileName(), dat, nio))
		//  {
		//  char *err = biffGetDone(DAT); // would be nice to free(err)
		//  itkExceptionMacro("Write: Error writing " 
		//                    << this->GetFileName() << ":\n" << err);
		//  }
		//
		//// Free the dat struct but don't touch dat->data
		//dat = datNix(dat);
		//nio = datIoStateNix(nio);
	}

} // end namespace itk
