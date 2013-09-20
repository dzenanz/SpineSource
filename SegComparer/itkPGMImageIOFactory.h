#ifndef __itkPGMImageIOFactory_h
#define __itkPGMImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class PGMImageIOFactory
 * \brief Create instances of PGMImageIO objects using an object factory.
 */
class ITK_EXPORT PGMImageIOFactory : public ObjectFactoryBase
{
public:  
  /** Standard class typedefs. */
  typedef PGMImageIOFactory       Self;
  typedef ObjectFactoryBase        Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Class methods used to interface with the registered factories. */
  virtual const char* GetITKSourceVersion(void) const;
  virtual const char* GetDescription(void) const;
  
  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PGMImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
    {
    PGMImageIOFactory::Pointer DatFactory = PGMImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(DatFactory);
    }

protected:
  PGMImageIOFactory();
  ~PGMImageIOFactory();

private:
  PGMImageIOFactory(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
  
  
} // end namespace itk

#endif
