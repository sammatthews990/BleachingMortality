
# For more info refer to https://github.com/DanTanAtAims/Reef-Monitoring/

using HTTP
using JSON
using DataFrames

# Base URL for the API
BASE_URL = "https://api.aims.gov.au/data-v2.0/10.25845/5c09bc4ff315c/"

# Endpoint to retrieve reef information
get_reef_info_url = joinpath(BASE_URL, "reef")

# Functions to construct URLs for different data types
get_photo_transect(name::String; domain_category::String = "reef") =
	joinpath(BASE_URL, "data?domain_name=$(HTTP.escapeuri(name))&domain_category=$(domain_category)&data_type=photo-transect")

get_manta_tow(name::String; domain_category::String = "reef") =
	joinpath(BASE_URL, "data?domain_name=$(HTTP.escapeuri(name))&domain_category=$(domain_category)&data_type=manta")

get_disturbances(name::String; aggregation::String = "reef") =
	joinpath(BASE_URL, "disturbance?reef=$(HTTP.escapeuri(name))&aggregation=$(aggregation)")

get_cots(name::String; domain_category::String = "reef") =
	joinpath(BASE_URL, "cots-by-domain?domain_category=$(domain_category)&domain_name=$(HTTP.escapeuri(name))")

# Function to retrieve and parse JSON data from a given URL
function fetch_json_data(url::String)
	response = HTTP.get(url)
	if response.status == 200
		return JSON.parse(String(response.body))
	else
		println("Failed to retrieve data from: $url")
		return nothing
	end
end

# Step 1: Retrieve the list of reefs
reef_info = fetch_json_data(get_reef_info_url)
reef_names = [reef["aims_reef_name"] for reef in reef_info]

# data = fetch_json_data(get_photo_transect("Agincourt Reef No.1"))
# df = DataFrame(data)
# @show typeof(data)
# @show typeof(data[1])
# @show data[1]
# data = fetch_json_data(get_disturbances("Agincourt Reef No.1"))
# df = DataFrame(data)
# @show typeof(data)
# @show typeof(data[1])
# @show data[1]

photo_transect_data = Dict{String, DataFrame}()
for reef_name in reef_names
	println("Fetching photo transect data for reef: $reef_name")
	url = get_photo_transect(reef_name)
	data = fetch_json_data(url)
	if data !== nothing && !isempty(data)
		photo_transect_data[reef_name] = DataFrame(data)
	end
end

combined_df = vcat([
	transform(df, :domain_name => (_ -> reef) => :reef_name)
	for (reef, df) in photo_transect_data
]...)

# Fetch and combine manta tow data
manta_data_dict = Dict{String, DataFrame}()
for reef in reef_names
	println("Fetching manta tow for: $reef")
	data = fetch_json_data(get_manta_tow(reef))
	if data !== nothing && !isempty(data)
		manta_data_dict[reef] = DataFrame(data)
	end
end

mantaHC_df = vcat([
	transform(df, :domain_name => (_ -> reef) => :reef_name)
	for (reef, df) in manta_data_dict
]...)

# Fetch and combine disturbance data
disturbance_data_dict = Dict{String, DataFrame}()
for reef in reef_names
	println("Fetching disturbances for: $reef")
	data = fetch_json_data(get_disturbances(reef))
	if data !== nothing && !isempty(data)
		disturbance_data_dict[reef] = DataFrame(data)
	end
end

# Get the union of all column names across all dataframes
all_cols = reduce(union, names.(values(disturbance_data_dict)))

# Add missing columns to each DataFrame before concatenation
disturbance_df = vcat([
	begin
		df_copy = copy(df)
		# Insert missing columns as `:colname => missing`
		for col in setdiff(all_cols, names(df_copy))
			insertcols!(df_copy, Symbol(col) => missing)
		end
		insertcols!(df_copy, :reef_name => reef)
	end
	for (reef, df) in disturbance_data_dict
]...)

# Fetch and combine COTS data
cots_data_dict = Dict{String, DataFrame}()
for reef in reef_names
	println("Fetching COTS data for: $reef")
	data = fetch_json_data(get_cots(reef))
	if data !== nothing && !isempty(data)
		cots_data_dict[reef] = DataFrame(data)
	end
end

cots_df = vcat([
	transform(df, :domain_name => (_ -> reef) => :reef_name)
	for (reef, df) in cots_data_dict
]...)

# Fetch photo transect data for each sector
sectors = unique([reef["a_sector"] for reef in reef_info])

sector_photo_data = Dict{String, DataFrame}()

for sector in sectors
	println("Fetching photo transect for sector: $sector")
	url = get_photo_transect(sector; domain_category = "sector")
	data = fetch_json_data(url)
	if data !== nothing && !isempty(data)
		sector_photo_data[sector] = DataFrame(data)
	end
end

sector_photo_df = vcat([
	insertcols!(copy(df), :sector_name => sector)
	for (sector, df) in sector_photo_data
]...)

# Fetch manta tow data for each sector
sector_manta_data = Dict{String, DataFrame}()

for sector in sectors
	println("Fetching manta tow for sector: $sector")
	url = get_manta_tow(sector; domain_category = "sector")
	data = fetch_json_data(url)
	if data !== nothing && !isempty(data)
		sector_manta_data[sector] = DataFrame(data)
	end
end

sector_manta_df = vcat([
	insertcols!(copy(df), :sector_name => sector)
	for (sector, df) in sector_manta_data
]...)

# Fetch cots manta tow data for each sector
sector_cots_data = Dict{String, DataFrame}()

for sector in sectors
	println("Fetching COTS data for sector: $sector")
	url = get_cots(sector; domain_category = "sector")
	data = fetch_json_data(url)
	if data !== nothing && !isempty(data)
		sector_cots_data[sector] = DataFrame(data)
	end
end

sector_cots_df = vcat([
	insertcols!(copy(df), :sector_name => sector)
	for (sector, df) in sector_cots_data
]...)
