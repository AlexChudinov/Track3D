#pragma once

#define WIN32_LEAN_AND_MEAN         
#define NOMINMAX
#include <windows.h>
#undef NOMINMAX

#ifndef _USE_LS_DLL_
#define LAPLACIAN_SOLVER_EXPORT __declspec(dllexport)
#else
#define LAPLACIAN_SOLVER_EXPORT __declspec(dllimport)
#endif // !_USE_LS_DLL_

EXTERN_C LAPLACIAN_SOLVER_EXPORT void testMsg(LPCWCH msg);

// TODO: Установите здесь ссылки на дополнительные заголовки, требующиеся для программы
