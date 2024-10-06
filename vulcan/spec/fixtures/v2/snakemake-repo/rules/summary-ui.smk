configfile: "config.yaml"

rule checker_ui:
    output:
        "output/ui_check.txt"
        # Since this is a UI job, we just specify the output, so snakemake can infer the dag

rule summary:
    input:
        poem_count="output/count_poem.txt",
        poem_count_2="output/count_poem_2.txt",
        arithmetic="output/arithmetic.txt",
        checker="output/check.txt",
        ui_checker="output/ui_check.txt"
    output:
        "output/summary.txt"
    run:
        with open(output[0], "w") as f:
            f.write(f"Count of poem.txt: {open(input.poem_count).read().strip()}\n")
            f.write(f"Count of poem_2.txt: {open(input.poem_count_2).read().strip()}\n")
            f.write(f"Result of arithmetic operation: {open(input.arithmetic).read().strip()}\n")
            f.write(f"Result of checker: {open(input.checker).read().strip()}\n")
            f.write(f"Result of UI checker: {open(input.ui_checker).read().strip()}\n")