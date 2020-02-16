package phoenixcenter.metaproteomics.entity;

import lombok.Data;
import lombok.NoArgsConstructor;

@Data
@NoArgsConstructor
public class UnipeptTaxon {
    private Integer id;
    private String name;
    private String rank;
}