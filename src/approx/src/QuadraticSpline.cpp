#include "QuadraticSpline.h"
#include <amp.h>  

#include <iostream>  

template <int TileSize>
floattype calcZ(floattype *tile_data, concurrency::tiled_index<TileSize> tidx) {
	
}

HRESULT CQuadraticSpline::Approximate(TApproximationParams *params) {
	size_t count;
	mEnumeratedLevels->GetLevelsCount(&count);
	TGlucoseLevel *levels;
	mEnumeratedLevels->GetLevels(&levels);	
	/*invalid levels*/
	if (count < 2) {
		return S_FALSE;
	}
	//size_t n = count - 1;
	floattype *z = (floattype*)malloc(count * sizeof(floattype));	
	//prepare z
	//t is datetime, y is level
	z[0] = 0;
	for (size_t i = 0; i < count - 1; i++) {
		//z{i+1} = -z{i} + 2*((y{i+1}-y{i})/(t{i+1}-t{i}))
		z[i + 1] = -z[i] + 2 * ((levels[i + 1].level - levels[i].level) / (levels[i + 1].datetime - levels[i].datetime));
		std::cout << z[i];
		//a{i} is y{i}, b{i} is z{i} and c{i}=(z{i+1}-z{i})/(2*(t{i+1}-t{i}))
		qParams.insert(std::pair<floattype,quadratic_params>(levels[i].datetime, quadratic_params(levels[i].level, z[i], (z[i+1]-z[i])/(2* (levels[i + 1].datetime - levels[i].datetime)))));
	}
	/*add last timestamp to check bounds in getLevels*/
	lastTime = levels[count - 1].datetime;
	//printf("aproximuju %d %f %f\n", count, cParams.begin()->first, lastTime);

#ifdef PARALLEL_AMP
	static const int TileSize = 64;
	std::vector<floattype> levelsV;
	std::vector<floattype> datetimesV;
	std::vector<quadratic_params> paramsV;
	floattype *zC = (floattype*)malloc(count * sizeof(floattype));
	//prepare vectors
	zC[0] = 0;
	for (size_t i = 0; i < count; i++) {
		levelsV.push_back(levels[i].level);
		datetimesV.push_back(levels[i].datetime);
		if (i == count - 1) {
			continue;
		}
		zC[i+1] = 2 * ((levels[i + 1].level - levels[i].level) / (levels[i + 1].datetime - levels[i].datetime));
		paramsV.push_back(quadratic_params(0, 0, 0));
	}
	// Create C++ AMP objects.  
	concurrency::array_view<floattype, 1> y(count, levelsV);
	concurrency::array_view<floattype, 1> x(count, datetimesV);
	concurrency::array_view<floattype, 1> zP(count, zC);
//	concurrency::array_view<quadratic_params, 1> pars(count, paramsV);
	//m.discard_data();
//	const auto compute_domain = pars.extent.tile<TileSize>().pad();
	const auto compute_domain = zP.extent.tile<TileSize>().pad();
	concurrency::array<floattype, 1> tile_sums(compute_domain / TileSize);
	concurrency::array_view<floattype, 1> tile_sums_vw(tile_sums);
	scan<64, scan_mode::inclusive>(concurrency::accelerator::get_auto_selection_view(),zP, zP, minus<floattype>());
	zP.synchronize();
	/*concurrency::parallel_for_each(
		// Define the compute domain, which is the set of threads that are created.  
		compute_domain,
		// Define the code to run on each thread on the accelerator.  
		[=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
		{
			const int gidx = tidx.global[0];
			const int lidx = tidx.local[0];
			//const int partial_data_length = pars.extent.size() - tidx.tile[0] * tidx.tile_extent[0];
			const int partial_data_length = zP.extent.size() - tidx.tile[0] * tidx.tile_extent[0];
			//tile static mem, z 0, x 1 y 2
			tile_static floattype tile_data[3][TileSize];
			tile_data[0][lidx] = 0;
			tile_data[1][lidx] = (lidx >= partial_data_length) ? 0 : x[gidx];
			tile_data[2][lidx] = (lidx >= partial_data_length) ? 0 : y[gidx];
			const floattype current_value = 0;// 2 * (tile_data[2][lidx] - tile_data[2][lidx - 1]) / (tile_data[1][lidx] - tile_data[1][lidx - 1]);
			tidx.barrier.wait_with_tile_static_memory_fence();

			//hnus
			for (int stride = 1; stride <= (TileSize / 2); stride *= 2)
			{
				if ((lidx + 1) % (stride * 2) == 0)
				{
					//z{i+1} = -z{i} + 2*((y{i+1}-y{i})/(t{i+1}-t{i}))
					tile_data[0][lidx] = -tile_data[0][lidx - stride] + 2 * (tile_data[2][lidx] - tile_data[2][lidx - 1]) / (tile_data[1][lidx] - tile_data[1][lidx - 1]);
				}
				tidx.barrier.wait_with_tile_static_memory_fence();
			}

			if (lidx == 0)
			{
				tile_data[0][TileSize - 1] = 0;
			}
			tidx.barrier.wait_with_tile_static_memory_fence();

			for (int stride = TileSize / 2; stride >= 1; stride /= 2)
			{
				if ((lidx + 1) % (stride * 2) == 0)
				{
					auto tmp = tile_data[0][lidx];
					tile_data[0][lidx] = tile_data[0][lidx] + tile_data[0][lidx - stride];
					tile_data[0][lidx - stride] = tmp;
				}
				tidx.barrier.wait_with_tile_static_memory_fence();
			}
			floattype ret = tile_data[0][TileSize - 1];
			//floattype z = calcZ<TileSize>(&(tile_data[0]), tidx);

			/*

			if (_Mode == scan_mode::inclusive)
			{
				tile_data[lidx] += current_value;
			}
			
			// This does not execute correctly on Warp accelerators for some reason. Steps Warp A & B do this instead.
			if (lidx == (TileSize - 1))
			{
				tile_sums_vw[tidx.tile[0]] = ret + current_value;
			}
			if (gidx < (zP.extent)[0])
			{
				zP[gidx] = tile_data[0][lidx];
			}
		});

	if (tile_sums_vw.extent[0] > TileSize)
	{
		scan<TileSize, amp_algorithms::scan_mode::exclusive>(accl_view, tile_sums_vw, tile_sums_vw, op);
	}
	else
	{
		concurrency::parallel_for_each(accl_view, compute_domain, [=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
		{
			const int gidx = tidx.global[0];
			const int lidx = tidx.local[0];
			const int partial_data_length = tile_sums_vw.extent.size() - tidx.tile[0] * tidx.tile_extent[0];

			tile_static T tile_data[TileSize];
			tile_data[lidx] = (lidx >= partial_data_length) ? 0 : tile_sums_vw[gidx];

			tidx.barrier.wait_with_tile_static_memory_fence();

			_details::scan_tile_exclusive<TileSize>(tile_data, tidx, amp_algorithms::plus<T>());

			tile_sums_vw[gidx] = tile_data[lidx];
			tidx.barrier.wait_with_tile_static_memory_fence();
		});
	}

	// 4. Add the tile results to the individual results for each tile.

	concurrency::parallel_for_each(accl_view, compute_domain, [=](concurrency::tiled_index<TileSize> tidx) restrict(amp)
	{
		const int gidx = tidx.global[0];

		if (gidx < output_view.extent[0])
		{
			output_view[gidx] += tile_sums_vw[tidx.tile[0]];
		}
	});
	*/

	
	zP.synchronize();
	for (int i = 0; i < count - 1; i++) {
		//printf("jedu %f %f\n", zP[i], z[i]);// mC[i]);
	}
#endif
	free(z);
	return S_OK;
}
HRESULT  CQuadraticSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder) {
	*filled = 0;
	/*Not aproximated yet*/
	if (qParams.empty()) {
		printf("err empty qparams\n");
		return S_FALSE;
	}
	//printf("vracim hodnoty %f %f %f", desiredtime, cParams.begin()->first, cParams.rbegin()->first);
	//cParams.end()->first
	//printf("test bond %f\n", cParams.lower_bound(desiredtime)->first);
	//printf("co se deje?\n");
	//std::map<floattype, cubic_params>::iterator itlow;
	//itlow = cParams.lower_bound(desiredtime);
	/*Invalid time (after end or too much before start)*/
	if (lastTime < desiredtime || desiredtime+(stepping*count) < qParams.begin()->first) {
		printf("err time\n");
		return S_FALSE;
	}
	//printf("prdelni testik %d\n", (cParams.find(desiredtime)));
	//std::map<floattype, cubic_params>::iterator it = cParams.find(desiredtime);
	/*Move time to start*/
	floattype currenttime = desiredtime;
	while (currenttime < qParams.begin()->first) {
		currenttime += stepping;
		count--;
	}
	/*Walk trough*/
	std::map<floattype, quadratic_params>::iterator it;
	floattype res;
	if (derivationorder == 0) {
		//printf("tu jsem1 %lu %f %f\n", *filled, currenttime, lastTime);
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);

			//workaround lower bound;
			if (it != qParams.begin() && it->first != currenttime)
			{
				it--;
			}
			/*workround for last value*/
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.a + it->second.b*(currenttime - it->first) + it->second.c*(currenttime - it->first)*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			//printf("tu jsem2 %lu\n", *filled);
		}
	}
	else if (derivationorder == 1){
		while (currenttime <= lastTime && count > 0) {
			it = qParams.lower_bound(currenttime);

			//workaround lower bound;
			if (it != qParams.begin() && it->first != currenttime)
			{
				it--;
			}
			//workround for last value
			if (it == qParams.end()) {
				it--;
			}
			res = it->second.b + 2*it->second.c*(currenttime - it->first);
			levels[(*filled)] = res;
			(*filled)++;
			currenttime += stepping;
			count--;
			//printf("tu jsem roblhl %lu\n", *filled);
		}
	}
	/*wrong derivative order*/
	else {
		return S_FALSE;
	}
	return S_OK;
}

CQuadraticSpline::CQuadraticSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {
}

CQuadraticSpline::~CQuadraticSpline() {
}
