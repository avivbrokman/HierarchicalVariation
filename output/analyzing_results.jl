using JSON

function fill_patch(patch, centers)
    return [centers[i] for i in patch]
end

function fill_partition(partition, centers)
    return [fill_patch(el, centers) for el in partition]
end

function make_aesthetic(quick_view)
    patch_strings = [join(sort(el), "/") for el in quick_view]
    return join(sort(patch_strings), ", ")
end

function fancy_print(file_path, p1)
    json_string = read(file_path, String)
    json_data = JSON.parse(json_string)
    centers = json_data["centers"]
    r_centers = round.(centers, digits = 3)
    partition = json_data["partition"]
    quick_view = fill_partition(partition, r_centers)
    print_view = make_aesthetic(quick_view)
    print("p1 = $p1: $print_view\n\n")
end


alpha1 = 0.6
p1s = 0.5:1/40:1
for p1 in p1s
    file_path = "output/finer1/hierarchical/fecundity4/delta0.1/bimodal/alpha1_$(alpha1)__p1_$p1/output.json"
    fancy_print(file_path, p1)
end

alpha1 = 0.9
p1s = 0.5:1/40:1
for p1 in p1s
    file_path = "output/finer2/hierarchical/fecundity4/delta0.1/bimodal/alpha1_$(alpha1)__p1_$p1/output.json"
    fancy_print(file_path, p1)
end

alpha1 = 0.9
p1s = 0.8:0.002:0.85
for p1 in p1s
    file_path = "output/finer1/hierarchical/fecundity4/delta0.1/bimodal/alpha1_$(alpha1)__p1_$p1/output.json"
    fancy_print(file_path, p1)
end
