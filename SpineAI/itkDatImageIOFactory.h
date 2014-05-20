#ifndef __itkDatImageIOFactory_h
#define __itkDatImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class DatImageIOFactory
 * \brief Create instances of DatImageIO objects using an object factory.
 */
class ITK_EXPORT DatImageIOFactory : public ObjectFactoryBase
{
public:  
  /** Standard class typedefs. */
  typedef DatImageIOFactory       Self;
  typedef ObjectFactoryBase        Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Class methods used to interface with the registered factories. */
  virtual const char* GetITKSourceVersion(void) const;
  virtual const char* GetDescription(void) const;
  
  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DatImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
    {
    DatImageIOFactory::Pointer DatFactory = DatImageIOFactory::New();
    ObjectFactoryBase::RegisterFactory(DatFactory);
    }

protected:
  DatImageIOFactory();
  ~DatImageIOFactory();

private:
  DatImageIOFactory(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
  
  
} // end namespace itk

#endif
