settings = Dict("barlen"=>20,
                "dt"=>0.1,
            )
progress_color = Dict("critical_line"=>:red,
                        "parallel_scans"=>:cyan,
                        "exploration"=>:yellow,
                )

set_progress_bar(N,name) = Progress(N, 
                                desc="Progress", 
                                color = progress_color[name],
                                barlen=settings["barlen"],
                                dt=settings["dt"]
                            )
