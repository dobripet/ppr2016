#pragma once

#include "CommonApprox.h"
#include "GlucoseLevels.h"
#include "..\..\common\rtl\referencedImpl.h"
#include <map>
#include <amp.h>

#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

typedef struct quadratic_params{
	floattype a;
	floattype b;
	floattype c;
	quadratic_params(floattype a, floattype b, floattype c) : a(a), b(b), c(c) {};
	quadratic_params() {};
} quadratic_params;


class CQuadraticSpline : public CCommonApprox, public virtual CReferenced {
protected:
	std::map<floattype, quadratic_params> qParams;
	floattype lastTime;
public:
	CQuadraticSpline(IGlucoseLevels *levels);
	virtual ~CQuadraticSpline();
	//dctor has to be virtual, even if it is empty, due to the inheritance by dominance	
	HRESULT Approximate(TApproximationParams * params);
	HRESULT GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype * levels, size_t * filled, size_t derivationorder);
};

#pragma warning( pop )


template <typename T>
class minus
{
public:
	T operator()(const T &a, const T &b) const restrict(cpu, amp)
	{
		return (a - b);
	}
};

template <typename T>
class plus
{
public:
	T operator()(const T &a, const T &b) const restrict(cpu, amp)
	{
		return (a + b);
	}
};

enum class scan_mode : int
{
	exclusive = 0,
	inclusive = 1
};

enum class scan_direction : int
{
	forward = 0,
	backward = 1
};


//----------------------------------------------------------------------------
// scan - C++ AMP implementation
//----------------------------------------------------------------------------
//
// References:
//
// "GPU Gems 3, chapter 39. Parallel Prefix Sum (Scan) with CUDA" http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html
//
// https://research.nvidia.com/sites/default/files/publications/nvr-2008-003.pdf
// https://sites.google.com/site/duanemerrill/ScanTR2.pdf
//
// TODO: There may be some better scan implementations that are described in the second reference. Investigate.
// TODO: Scan only supports Rank of 1.
// TODO: Scan does not support forwards/backwards.

static const int scan_default_tile_size = 512;

template <int TileSize, typename _BinaryOp, typename T>
inline T scan_tile_exclusive(T* const tile_data, concurrency::tiled_index<TileSize> tidx, const _BinaryOp& op, const int partial_data_length) restrict(amp)
{
	const int lidx = tidx.local[0];

	for (int stride = 1; stride <= (TileSize / 2); stride *= 2)
	{
		if ((lidx + 1) % (stride * 2) == 0)
		{
			tile_data[lidx] = op(tile_data[lidx], tile_data[lidx - stride]);
		}
		tidx.barrier.wait_with_tile_static_memory_fence();
	}

	if (lidx == 0)
	{
		tile_data[TileSize - 1] = 0;
	}
	tidx.barrier.wait_with_tile_static_memory_fence();

	for (int stride = TileSize / 2; stride >= 1; stride /= 2)
	{
		if ((lidx + 1) % (stride * 2) == 0)
		{
			auto tmp = tile_data[lidx];
			tile_data[lidx] = op(tile_data[lidx], tile_data[lidx - stride]);
			tile_data[lidx - stride] = tmp;
		}
		tidx.barrier.wait_with_tile_static_memory_fence();
	}
	return tile_data[TileSize - 1];
}
//helpers
template <int N, typename InputIndexableView>
inline int tile_partial_data_size(const InputIndexableView& arr, concurrency::tiled_index<N> tidx) restrict(amp)
{
	return arr.extent.size() - tidx.tile[0] * tidx.tile_extent[0];
}
template <typename InputIndexableView, int N>
inline void padded_write(InputIndexableView& arr, const concurrency::index<N> idx, const typename InputIndexableView::value_type &value) restrict(cpu, amp)
{
	if (arr.extent.contains(idx))
	{
		arr[idx] = value;
	}
}
template <typename InputIndexableView>
inline void padded_write(InputIndexableView& arr, const int idx, const typename InputIndexableView::value_type &value) restrict(cpu, amp)
{
	padded_write<InputIndexableView, 1>(arr, concurrency::index<1>(idx), value);
}
template <int TileSize, scan_mode _Mode, typename _BinaryFunc, typename InputIndexableView>
inline void scan(const concurrency::accelerator_view& accl_view, const InputIndexableView& input_view, InputIndexableView& output_view, const _BinaryFunc& op)
{
	typedef InputIndexableView::value_type T;

	const auto compute_domain = output_view.extent.tile<TileSize>().pad();
	concurrency::array<T, 1> tile_results(compute_domain / TileSize, accl_view);
	concurrency::array_view<T, 1> tile_results_vw(tile_results);

	// 1 & 2. Scan all tiles and store results in tile_results.

	concurrency::parallel_for_each(accl_view, compute_domain, [=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
	{
		const int gidx = tidx.global[0];
		const int lidx = tidx.local[0];
		const int partial_data_length = tile_partial_data_size(output_view, tidx);

		tile_static T tile_data[TileSize];
		tile_data[lidx] = (lidx >= partial_data_length) ? 0 : input_view[gidx];
		const T current_value = tile_data[lidx];
		tidx.barrier.wait_with_tile_static_memory_fence();

		auto val = scan_tile_exclusive<TileSize>(tile_data, tidx, plus<T>(), partial_data_length);
		if (_Mode == scan_mode::inclusive)
		{
			tile_data[lidx] += current_value;
		}

		if (lidx == (TileSize - 1))
		{
			tile_results_vw[tidx.tile[0]] = val + current_value;
		}
		padded_write(output_view, gidx, tile_data[lidx]);
	});

	// 3. Scan tile results.

	if (tile_results_vw.extent[0] > TileSize)
	{
		scan<TileSize, scan_mode::exclusive>(accl_view, tile_results_vw, tile_results_vw, op);
	}
	else
	{
		concurrency::parallel_for_each(accl_view, compute_domain, [=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
		{
			const int gidx = tidx.global[0];
			const int lidx = tidx.local[0];
			const int partial_data_length = tile_partial_data_size(tile_results_vw, tidx);

			tile_static T tile_data[TileSize];
			tile_data[lidx] = tile_results_vw[gidx];
			tidx.barrier.wait_with_tile_static_memory_fence();

			scan_tile_exclusive<TileSize>(tile_data, tidx, plus<T>(), partial_data_length);

			tile_results_vw[gidx] = tile_data[lidx];
			tidx.barrier.wait_with_tile_static_memory_fence();
		});
	}

	// 4. Add the tile results to the individual results for each tile.

	concurrency::parallel_for_each(accl_view, compute_domain, [=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
	{
		const int gidx = tidx.global[0];

		if (gidx < output_view.extent[0])
		{
			output_view[gidx] += tile_results_vw[tidx.tile[0]];
		}
	});
}


