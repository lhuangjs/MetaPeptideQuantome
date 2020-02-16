package phoenixcenter.metaproteomics.entity;

import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.List;

@Data
@NoArgsConstructor
public class UnipeptLCA {
    private String sequence;
    private Integer lcaId;
    private String lcaName;
    private String lcaRank;
    private List<Integer> lineageIds;
    private List<String> lineageNames;
}
