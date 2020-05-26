#pragma once

using namespace System;
namespace FMatVecClr {

    template<class T>
    public ref class ManagedObject
    {
    protected:
        T* m_Instance;
    public:
        ManagedObject(T* instance)
            : m_Instance(instance)
        {
        }
        virtual ~ManagedObject()
        {
            if (m_Instance != nullptr)
            {
                delete m_Instance;
            }
        }
        !ManagedObject()
        {
            if (m_Instance != nullptr)
            {
                delete m_Instance;
            }
        }
        T* GetInstance()
        {
            return m_Instance;
        }
    };
}

//using namespace System::Runtime::InteropServices;
//static const char* string_to_char_array(String^ string)
//{
//    const char* str = (const char*)(Marshal::StringToHGlobalAnsi(string)).ToPointer();
//    return str;
//}