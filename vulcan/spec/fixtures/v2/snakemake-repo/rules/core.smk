configfile: "config.yaml"

rule count:
    input:
        # Any jobs that use a config as input and are file paths, must make sure the config values
        # match the vulcan_config.yaml
        poem1=config["poem"],
        poem2=config["poem_2"]
    output:
        poem_count="output/count_poem.txt",
        poem_count_2="output/count_poem_2.txt"
    params:
        count_bytes=config["count_bytes"],
        count_chars=config["count_chars"]
    shell:
        """
        if {params.count_bytes}; then
            wc_option="-c"
        elif {params.count_chars}; then
            wc_option="-m"
        else
            wc_option="-w"
        fi
        wc $wc_option {input.poem1} | awk '{{print $1}}' > {output.poem_count}
        wc $wc_option {input.poem2} | awk '{{print $1}}' > {output.poem_count_2}
        """

rule arithmetic:
    input:
        poem_count="output/count_poem.txt",
        poem_count_2="output/count_poem_2.txt"
    output:
        "output/arithmetic.txt"
    params:
        add=config["add"],
        add_and_multiply_by=config["add_and_multiply_by"]
    run:
        with open(input.poem_count, "r") as f:
            count1 = int(f.read().strip())
        with open(input.poem_count_2, "r") as f:
            count2 = int(f.read().strip())

        if params.add:
            result = count1 + count2
        else:
            result = (count1 + count2) * params.add_and_multiply_by

        with open(output[0], "w") as f:
            f.write(str(result))

rule checker:
    input:
        "output/arithmetic.txt"
    output:
        "output/check.txt"
    run:
        with open(input[0], "r") as f:
            if not f.read().strip():
                raise ValueError("No value found in arithmetic.txt")
        with open(output[0], "w") as f:
            f.write("Value exists in arithmetic.txt\n")